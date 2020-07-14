import uproot
import numpy as np
import numba
import awkward1 as ak
import os
import optparse

from data import getData
from plot_utils import plotHist, plotCollection, plotEfficiency, plotROC, plotGraphs, combinePDFs
from truth_utils import truth_link, dr_match
from plot_config import reco_pl,truth_pl

def analyze(opts, args):
    
    cms_events = getData(opts)
    nEvents=len(cms_events)
            
    # histogram all attributes for input quantities
    if opts.drawInputHistograms:
        for collection in cms_events.columns:
            plotCollection(cms_events[collection], collection, xtitle=collection, outDir=opts.odir+"/diagnostic/input_collections")    

    # common ficucial region selection
    pt_truth_lo = 2.5
    pt_truth_hi = 4.5
    pt_reco_lo = 2
    pt_reco_hi = 5
            
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
    # fiducial region selection
    truth_electrons = truth_electrons[ (truth_electrons.pt > pt_truth_lo) & (truth_electrons.pt < pt_truth_hi) ]
    
    if opts.drawTruthElectrons:
        plotCollection(truth_electrons, "genElectrons", xtitle="genElectrons", outDir=opts.odir+"/diagnostic/all_truth_eles")    
    
    # DEFINE THE RECO ELECTRONS
    #reco_electrons = cms_events['electrons']
    reco_electrons = cms_events['softElectrons']
    reco_cuts = {}
    reco_cuts["all"] = (reco_electrons.pt > pt_reco_lo) & (reco_electrons.pt < pt_reco_hi) 
   # reco_cuts["presel"] = (reco_cuts["all"] & 
   #                        (np.abs(reco_electrons.dxy)< 0.1 ) & (np.abs(reco_electrons.dz) < 15 ) & 
   #                        (reco_electrons.ip3d < 5 ) & (reco_electrons.trkRelIso < 2 ) &
   #                        (reco_electrons.mvaId>-1) & (reco_electrons.ptBiased>-1)
   #                       )
   # reco_cuts["doptimization.dxy"] = (reco_cuts["presel"] & (np.abs(reco_electrons.dxy) < 0.025 ))
   # reco_cuts["doptimization.dz"] = (reco_cuts["presel"] & (np.abs(reco_electrons.dz) < 0.1 ))
   # reco_cuts["doptimization.ip3d"] = (reco_cuts["presel"] & (reco_electrons.ip3d < 2 ))
   # reco_cuts["doptimization.mvaId"] = (reco_cuts["presel"] & (reco_electrons.mvaId > 2.5 ))
   # reco_cuts["doptimization.fBrem"] = (reco_cuts["presel"] & (reco_electrons.fBrem > 0.05 ))
   # reco_cuts["doptimization.ptBiased"] = (reco_cuts["presel"] & (reco_electrons.ptBiased > 2.5 ))
   # reco_cuts["doptimization.unBiaed"] = (reco_cuts["presel"] & (reco_electrons.unBiased > 2.5 ))
   # reco_cuts["doptimization.trkRelIso"] = (reco_cuts["presel"] & (reco_electrons.trkRelIso < 0.3 ))
   # reco_cuts["doptimization.sip3d"] = (reco_cuts["presel"] & (reco_electrons.sip3d < 2.5 ))

    reco_cuts["cut1_mvaId"] = (reco_cuts["all"] & (reco_electrons.mvaId > 3. ))
    reco_cuts["cut2_dxy"] = (reco_cuts["cut1_mvaId"] & (reco_electrons.dxy < 0.02) & (reco_electrons.dxy > -0.02))
    reco_cuts["cut3_dz"] = (reco_cuts["cut2_dxy"] & (reco_electrons.dz < 8) & (reco_electrons.dxy > -8))
    reco_cuts["cut4_trkRelIso"] = (reco_cuts["cut3_dz"] & (reco_electrons.trkRelIso<1))

    # Direct which ROCs to produce for each variable
    # also must describe what values are 'signal-like' (hi, low,
    # low absolute value [abslo], or high absolute values [abshi] )
    roc_config = {
        "all": {
            'dxy'         : "abslo",
            'dz'          : "abslo",
            'ip3d'        : "abslo",
            'sip3d'       : "lo",
            'trkRelIso'   : "lo",
            'mvaId'       : "hi",
            'ptBiased'    : "hi",
            'unBiased'    : "hi",
            },
        "cut1_mvaId": {
            'dxy'         : "abslo",
            'dz'          : "abslo",
            'ip3d'        : "abslo",
            'sip3d'       : "lo",
            'trkRelIso'   : "lo",
            'mvaId'       : "hi",
            'ptBiased'    : "hi",
            'unBiased'    : "hi",
            },
        "cut2_dxy": {
            'dxy'         : "abslo",
            'dz'          : "abslo",
            'ip3d'        : "abslo",
            'sip3d'       : "lo",
            'trkRelIso'   : "lo",
            'mvaId'       : "hi",
            'ptBiased'    : "hi",
            'unBiased'    : "hi",
            },
        "cut3_dz": {
            'dxy'         : "abslo",
            'dz'          : "abslo",
            'ip3d'        : "abslo",
            'sip3d'       : "lo",
            'trkRelIso'   : "lo",
            'mvaId'       : "hi",
            'ptBiased'    : "hi",
            'unBiased'    : "hi",
            },
        "cut4_trkRelIso": {
            'dxy'         : "abslo",
            'dz'          : "abslo",
            'ip3d'        : "abslo",
            'sip3d'       : "lo",
            'trkRelIso'   : "lo",
            'mvaId'       : "hi",
            'ptBiased'    : "hi",
            'unBiased'    : "hi",
            },
        # "presel": {
        #     'dxy'         : "abslo",
        #     'dz'          : "abslo",
        #     'ip3d'        : "abslo",
        #     'sip3d'       : "lo",
        #     'trkRelIso'   : "lo",
        #     'mvaId'       : "hi",
        #     'ptBiased'    : "hi",
        #     'unBiased'    : "hi",
        #     },
    }

    if opts.drawRecoElectrons:
        plotCollection(reco_electrons, "recoElectrons", xtitle="recoElectrons", outDir=opts.odir+"/diagnostic/all_reco_eles")    
    
    # BEGIN THE ANALYSIS
    evt_mask = ak.num(truth_electrons)==2
    
    efficiencies={}
    efficiencies_lo={}
    nFakes={}
    nFakes_lo={}
    for cut_name in reco_cuts:
        reco_mask = reco_cuts[cut_name] #(reco_electrons)
        # reco_mask = reco_cuts[cut_name]
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
            normalizeVars = False
            if cut_name == "all": normalizeVars = True
            plotCollection([matched_reco,unmatched_reco], cut_name+"_matching_reco", leg=["matched","unmatched"], plotlist=reco_pl,
                           outDir=opts.odir+"/diagnostic/match_comparison/reco_"+cut_name, normAttrs=normalizeVars, profile=True)
            if cut_name in roc_config:
                maxeff = ak.count(truth_electrons)
                if maxeff: maxeff = ak.count(matched_truth) / maxeff
                rocs={}
                for var in roc_config[cut_name]:
                    signal = ak.to_numpy(ak.flatten(matched_reco[var]))
                    background = ak.to_numpy(ak.flatten(unmatched_reco[var]))
                    roc = plotROC(var, signal, background, maxeff=maxeff, nEvts=nEvents, signalLike=roc_config[cut_name][var],
                            outDir=opts.odir+"/diagnostic/match_comparison/reco_"+cut_name+"/roc")
                    rocs[var]=roc
                plotGraphs("overlay_1d", [rocs[v] for v in rocs], leg=[v for v in rocs],
                           xtitle="Efficiency", ytitle="Fake multiplicity per evt",
                           outDir=opts.odir+"/diagnostic/match_comparison/reco_"+cut_name+"/roc")
                # passVals = ak.to_numpy(ak.flatten(matched_truth.pt))
                # totVals = ak.to_numpy(ak.flatten(truth_electrons.pt))
                # roc_cfg=roc_config[cut_name]
                # normalizeVars = False
        
        # display matching efficiencies
        passVals = ak.to_numpy(ak.flatten(matched_truth.pt))
        totVals = ak.to_numpy(ak.flatten(truth_electrons.pt))
        eff = plotEfficiency("eff_pt_"+cut_name, passVals=passVals, totVals=totVals, lims=(0,10), nbins=20,
                             xtitle="truth electron p_T [GeV]", outDir=opts.odir+"/efficiencies")[0]
        eff_lo = plotEfficiency("eff_pt_lo_"+cut_name, passVals=passVals, totVals=totVals, lims=(0.6,3), nbins=12,
                                xtitle="truth electron p_T [GeV]", outDir=opts.odir+"/efficiencies")[0]
        efficiencies[cut_name]=eff
        efficiencies_lo[cut_name]=eff_lo
    
        # record number of signal electrons
        nfakes = ak.num(signal_electrons) - ak.num(matched_reco)

        fakes = plotHist("n_signalElectrons_"+cut_name, vals=nfakes, lims=(-0.5,159.5), nbins=80,
                         xtitle="signal electron multiplicity", outDir=opts.odir+"/nFakes")

        fakes_lo = plotHist("n_signalElectrons_lo_"+cut_name, vals=nfakes, lims=(-0.5,19.5), nbins=20,
                            xtitle="signal electron multiplicity", outDir=opts.odir+"/nFakes")

        nFakes[cut_name]=fakes
        nFakes_lo[cut_name]=fakes_lo
    
    plotEfficiency("eff_pt",
                   effs=[efficiencies[cut] for cut in reco_cuts],
                   leg=[cut for cut in reco_cuts],
                   xtitle="truth electron p_T [GeV]",
                   outDir=opts.odir+"/final_comparisons")
    plotEfficiency("effLo_pt",
                   effs=[efficiencies_lo[cut] for cut in reco_cuts],
                   leg=[cut for cut in reco_cuts],
                   xtitle="truth electron p_T [GeV]",
                   outDir=opts.odir+"/final_comparisons")
    plotHist("nfakes",
             hists=[nFakes[cut] for cut in reco_cuts],
             leg=[cut for cut in reco_cuts], showMean=True,
             xtitle="fake multiplicity",
             outDir=opts.odir+"/final_comparisons")
    plotHist("nfakesLo",
             hists=[nFakes_lo[cut] for cut in reco_cuts],
             leg=[cut for cut in reco_cuts], showMean=True,
             xtitle="fake multiplicity",
             outDir=opts.odir+"/final_comparisons")
        
    combinePDFs()
    
    
    
    
if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-i',"--input", type="string", default = '', help="path to input file")
    parser.add_option('-o',"--odir", type="string", default = 'output_plots/', help="plots directory")
    parser.add_option('-b',"--bigInput", action='store_true', default = False, help="use large input file")
    parser.add_option("--drawInputHistograms", action='store_true', default = False, help="histogram ALL input collections")
    parser.add_option("--drawTruthElectrons", action='store_true', default = False, help="histogram selected truth electrons")
    parser.add_option("--drawRecoElectrons", action='store_true', default = False, help="histogram selected reco electrons")
    parser.add_option("--drawMatchedCollections", action='store_true', default = False, help="histogram all matched/unmatched truth and reco eles")
    parser.add_option("--drawMatchedReco", action='store_true', default = False, help="histogram all matched/unmatched reco eles")
    parser.add_option("--drawMatchedTruth", action='store_true', default = False, help="histogram all matched/unmatched truth eles")
    (options, args) = parser.parse_args()
    analyze(options, args)
