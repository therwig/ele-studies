import uproot
import numpy as np
import numba
import awkward1 as ak
import os
import optparse

from data import getData
from plot_utils import plotHist, plotCollection, plotEfficiency, combinePDFs
from truth_utils import truth_link, dr_match

def analyze(opts, args):
    
    cms_events = getData(opts)
            
    # histogram all attributes for input quantities
    if opts.drawInputHistograms:
        for collection in cms_events.columns:
            plotCollection(cms_events[collection], collection, xtitle=collection, outDir=opts.odir)    
    
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
        plotCollection(truth_electrons, "genElectrons", xtitle="genElectrons", outDir=opts.odir)    
    
    # DEFINE THE RECO ELECTRONS        
    #reco_electrons = cms_events['electrons']
    reco_electrons = cms_events['softElectrons']
    reco_cuts = {
        "all": reco_electrons.pt > -1,
        "looseMVA": reco_electrons.mvaId > 0,
        "tightMVA": reco_electrons.mvaId > 5,
        }
    
    if opts.drawRecoElectrons:
        plotCollection(reco_electrons, "recoElectrons", xtitle="recoElectrons", outDir=opts.odir)    
    
    # BEGIN THE ANALYSIS
    evt_mask = ak.num(truth_electrons)==2

    efficiencies={}
    nFakes={}
    for cut_name in reco_cuts:
        reco_mask = reco_cuts[cut_name]
        signal_electrons = reco_electrons[reco_mask]
        
        match_builder = dr_match(truth_electrons, signal_electrons, ak.ArrayBuilder())
        match_extension = match_builder.snapshot()
        matching_truth_mask = match_extension.truth_to_reco_index >= 0
        matching_reco_mask = match_extension.truth_to_reco_index[matching_truth_mask]
        #matching_reco_mask = match_extension.reco_to_truth_index >= 0
        
        matched_truth = truth_electrons[matching_truth_mask]
        unmatched_truth = truth_electrons[~matching_truth_mask]
        matched_reco = signal_electrons[matching_reco_mask]
        unmatched_reco = signal_electrons[~matching_reco_mask]
        
        if opts.drawMatchedCollections:
            plotCollection(matched_truth,   "matched_truth_ele",   xtitle="matched truth electron"  ,outDir=opts.odir)
            plotCollection(unmatched_truth, "unmatched_truth_ele", xtitle="unmatched truth electron",outDir=opts.odir)
            plotCollection(matched_reco,    "matched_truth_reco",  xtitle="matched reco electron"   ,outDir=opts.odir)
            plotCollection(unmatched_reco,  "unmatched_truth_reco",xtitle="unmatched reco electron" ,outDir=opts.odir)
        if opts.drawMatchedTruth:
            plotCollection([matched_truth,unmatched_truth], cut_name+"_matching_truth", leg=["matched","unmatched"],outDir=opts.odir)
        if opts.drawMatchedReco:
            plotCollection([matched_reco,unmatched_reco], cut_name+"_matching_reco", leg=["matched","unmatched"],outDir=opts.odir)
        
        # display matching efficiencies
        passVals = ak.to_list(ak.flatten(matched_truth.pt))
        totVals = ak.to_list(ak.flatten(truth_electrons.pt))
        eff = plotEfficiency("eff_pt_"+cut_name, passVals=passVals, totVals=totVals, lims=(0,20), nbins=20, xtitle="truth electron p_T [GeV]", outDir=opts.odir)[0]
        efficiencies[cut_name]=eff

        # record number of signal electrons
        fakes = plotHist("n_signalElectrons_"+cut_name, vals=ak.num(signal_electrons), xtitle="signal electron multiplicity", outDir=opts.odir, lims=(-0.5,99.5), nbins=50)
        nFakes[cut_name]=fakes

    plotEfficiency("eff_pt",
                   effs=[efficiencies[cut] for cut in reco_cuts],
                   leg=[cut for cut in reco_cuts],
                   xtitle="truth electron p_T [GeV]",
                   outDir=opts.odir)
    plotHist("nfakes",
             hists=[nFakes[cut] for cut in reco_cuts],
             leg=[cut for cut in reco_cuts],
             xtitle="fake multiplicity",
             outDir=opts.odir)
        
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
    
