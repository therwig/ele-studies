import uproot
import numpy as np
import numba
import awkward1 as ak
import os
import optparse

from data import getData
from plot_utils import plotHist, plotCollection, plotEfficiency, combinePDFs
from truth_utils import truth_link, dr_match
from plot_config import reco_pl,truth_pl

def analyze(opts, args):
    
    cms_events = getData(opts)
            
    # histogram all attributes for input quantities
    if opts.drawInputHistograms:
        for collection in cms_events.columns:
            plotCollection(cms_events[collection], collection, xtitle=collection, outDir=opts.odir+"/diagnostic/input_collections")    
    
    # DEFINE THE TRUTH ELECTRONS        
    # derived array with extra truth information
    truth_builder = truth_link(cms_events['genParticles'], ak.ArrayBuilder())
    truth_extension = truth_builder.snapshot()
    ele_mask = np.abs(cms_events['genParticles']['pdgId']) == 11
    first_mask = truth_extension['isFirst']
    last_mask = np.abs(cms_events['genParticles']['status']) == 1 # equivalent to last
    z_mask = truth_extension['motherPdgId'] == 23
    n2_mask = truth_extension['motherPdgId'] == 1000023
    gen_ele_mask = ele_mask & (z_mask | n2_mask) & last_mask
    truth_electrons = cms_events['genParticles'][gen_ele_mask]
    
    if opts.drawTruthElectrons:
        plotCollection(truth_electrons, "genElectrons", xtitle="genElectrons", outDir=opts.odir+"/diagnostic/all_truth_eles")    
    
    # DEFINE THE RECO ELECTRONS        
    #reco_electrons = cms_events['electrons']
    reco_electrons = cms_events['softElectrons']
    reco_cuts = {
        "all": reco_electrons.pt > -1,
        "looseMVA": reco_electrons.mvaId > 0,
        "tightMVA": reco_electrons.mvaId > 5,
        "dominic_special_ip3d"   :((reco_electrons.sip3d       < 1000)   &
                                   (reco_electrons.ip3d        < 8)      ),
        "dominic_special_dxyz"   :((np.abs(reco_electrons.dz)  < 0.1)    & 
                                   (reco_electrons.dzErr       < 0.5)    & 
                                   (np.abs(reco_electrons.dxy) < 0.05)   & 
                                   (reco_electrons.dxyErr      < 0.5)    ),
 #      "dominic_special_fBrem"  :  reco_electrons.fBrem       > - ,
        "dominic_special_Flav"   :  reco_electrons.GenPartFlav < 0.5     , 
        "dominic_special_mvaId"  :((reco_electrons.mvaId       > 2.5)    &
                                   (reco_electrons.mvaId       < 20)     ),
 #      "dominic_special_pt"     :  reco_electrons.pt          > - ,
        "dominic_special_ptBias" :((reco_electrons.ptBiased    > 2.5)    &
                                   (reco_electrons.ptBiased    < 20)     ),
 #      "dominic_special_trk"    :  reco_electrons.trkRelIso   < - ,
        "dominic_special_unBias" :((reco_electrons.unBiased    > 0)      &
                                   (reco_electrons.unBiased    < 5)      ), 
    }
    
    if opts.drawRecoElectrons:
        plotCollection(reco_electrons, "recoElectrons", xtitle="recoElectrons", outDir=opts.odir+"/diagnostic/all_reco_eles")    
    
    # BEGIN THE ANALYSIS
    evt_mask = ak.num(truth_electrons)==2
    
    efficiencies={}
    nFakes={}
    for cut_name in reco_cuts:
        reco_mask = reco_cuts[cut_name]
        signal_electrons = reco_electrons[reco_mask]

        debugMatching=False
        match_builder = dr_match(truth_electrons, signal_electrons, ak.ArrayBuilder(), doReco=False, debug=debugMatching)
        match_extension = match_builder.snapshot()
        matching_truth_mask = match_extension.truth_to_reco_index >= 0
        # matching_reco_mask = match_extension.truth_to_reco_index[matching_truth_mask]
    
        match_builderR = dr_match(truth_electrons, signal_electrons, ak.ArrayBuilder(), doReco=True, debug=debugMatching)
        match_extensionR = match_builderR.snapshot()
        matching_reco_mask = match_extensionR.reco_to_truth_index >= 0
        nonmatching_reco_mask = match_extensionR.reco_to_truth_index < 0
        
        
        matched_truth = truth_electrons[matching_truth_mask]
        unmatched_truth = truth_electrons[~matching_truth_mask]
        matched_reco = signal_electrons[matching_reco_mask]
        unmatched_reco = signal_electrons[~matching_reco_mask]
        # unmatched_reco = signal_electrons[nonmatching_reco_mask]
        
        if opts.drawMatchedCollections:
            plotCollection(matched_truth,   "matched_truth_ele",   xtitle="matched truth electron"  ,outDir=opts.odir+"/diagnostic/matched_collections/"+cut_name)
            plotCollection(unmatched_truth, "unmatched_truth_ele", xtitle="unmatched truth electron",outDir=opts.odir+"/diagnostic/matched_collections/"+cut_name)
            plotCollection(matched_reco,    "matched_reco_ele",  xtitle="matched reco electron"   ,outDir=opts.odir+"/diagnostic/matched_collections/"+cut_name)
            plotCollection(unmatched_reco,  "unmatched_reco_ele",xtitle="unmatched reco electron" ,outDir=opts.odir+"/diagnostic/matched_collections/"+cut_name)
        if opts.drawMatchedTruth:
            plotCollection([matched_truth,unmatched_truth], cut_name+"_matching_truth", leg=["matched","unmatched"], plotlist=truth_pl,
                           outDir=opts.odir+"/diagnostic/match_comparison/truth_"+cut_name, normAttrs=True)
        if opts.drawMatchedReco:
            plotCollection([matched_reco,unmatched_reco], cut_name+"_matching_reco", leg=["matched","unmatched"], plotlist=reco_pl,
                           outDir=opts.odir+"/diagnostic/match_comparison/reco_"+cut_name, normAttrs=True, profile=True)
        
        # display matching efficiencies
        passVals = ak.to_list(ak.flatten(matched_truth.pt))
        totVals = ak.to_list(ak.flatten(truth_electrons.pt))
        eff = plotEfficiency("eff_pt_"+cut_name, passVals=passVals, totVals=totVals, lims=(0,20), nbins=20, xtitle="truth electron p_T [GeV]", outDir=opts.odir+"/efficiencies")[0]
        efficiencies[cut_name]=eff
    
        # record number of signal electrons
        nfakes = ak.num(signal_electrons) - ak.num(matched_reco)
        fakes = plotHist("n_signalElectrons_"+cut_name, vals=nfakes, xtitle="signal electron multiplicity", outDir=opts.odir+"/nFakes", lims=(-0.5,159.5), nbins=80)
        nFakes[cut_name]=fakes
    
    plotEfficiency("eff_pt",
                   effs=[efficiencies[cut] for cut in reco_cuts],
                   leg=[cut for cut in reco_cuts],
                   xtitle="truth electron p_T [GeV]",
                   outDir=opts.odir+"/final_comparisons")
    plotHist("nfakes",
             hists=[nFakes[cut] for cut in reco_cuts],
             leg=[cut for cut in reco_cuts],
             xtitle="fake multiplicity",
             outDir=opts.odir+"/final_comparisons")
        
    combinePDFs()
    
    
    
    
if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-i',"--input", type="string", default = '', help="path to input file")
    parser.add_option('-o',"--odir", type="string", default = 'output_plots/', help="plots directory")
    parser.add_option("--drawInputHistograms", action='store_true', default = False, help="histogram ALL input collections")
    parser.add_option("--drawTruthElectrons", action='store_true', default = False, help="histogram selected truth electrons")
    parser.add_option("--drawRecoElectrons", action='store_true', default = False, help="histogram selected reco electrons")
    parser.add_option("--drawMatchedCollections", action='store_true', default = False, help="histogram all matched/unmatched truth and reco eles")
    parser.add_option("--drawMatchedReco", action='store_true', default = False, help="histogram all matched/unmatched reco eles")
    parser.add_option("--drawMatchedTruth", action='store_true', default = False, help="histogram all matched/unmatched truth eles")
    (options, args) = parser.parse_args()
    analyze(options, args)
