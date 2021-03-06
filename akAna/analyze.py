import uproot
import numpy as np
import numba
import awkward1 as ak
import os
import optparse

from data import getData
from plot_utils import plotHist, plotCollection, plotEfficiency, plotROC, plotGraphs, combinePDFs
from truth_utils import truth_link, dr_match, truth_reco_match
from plot_config import reco_pl,truth_pl

# def absLess(x,cut):
#     return (x < cut) & (-x > cut)

def analyze(opts, args):
    
    cms_events = getData(opts)
    nEvents=len(cms_events)
            
    # histogram all attributes for input quantities
    if opts.drawInputHistograms:
        for collection in cms_events.columns:
            plotCollection(cms_events[collection], collection, xtitle=collection, outDir=opts.odir+"/diagnostic/input_collections")    

    # common ficucial region selection
    pt_truth_lo = 0
    pt_truth_hi = 15 #4.5
    pt_reco_lo = 0 #
    pt_reco_hi = 15 #

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
    tight_IP = (np.abs(reco_electrons.dz) < 0.2) & (reco_electrons.sip3d < 3)
    reco_cuts["ip"] = reco_cuts["all"] & tight_IP
    signal_electrons = reco_electrons[ reco_cuts["ip"] ]
    
    # Direct which ROCs to produce for each variable
    # also must describe what values are 'signal-like' (hi, low,
    # low absolute value [abslo], or high absolute values [abshi] )
    roc_config = {
        'dxy'         : "abslo",
        'dz'          : "abslo",
        'ip3d'        : "abslo",
        'sip3d'       : "lo",
        'trkRelIso'   : "lo",
        'mvaId'       : "hi",
        'ptBiased'    : "hi",
        'unBiased'    : "hi",
    }

    if opts.drawRecoElectrons:
        plotCollection(reco_electrons, "recoElectrons", xtitle="recoElectrons", outDir=opts.odir+"/diagnostic/all_reco_eles")    
    
    # BEGIN THE ANALYSIS
    met = cms_events['met']
    metPhi = cms_events['metPhi']
    mN1 = cms_events['genParticles'][ cms_events['genParticles'].pdgId==1000022 ][:,0].mass
    mN2 = cms_events['genParticles'][ cms_events['genParticles'].pdgId==1000023 ][:,0].mass
    dM = mN2 - mN1
    plotHist("mN1", vals=mN1, xtitle="mN1", outDir=opts.odir+"/truthAna")
    plotHist("mN2", vals=mN2, xtitle="mN2", outDir=opts.odir+"/truthAna")
    plotHist("dM", vals=dM, xtitle="dM", outDir=opts.odir+"/truthAna")
    plotHist("nEle", vals=ak.num(truth_electrons), xtitle="n truth electrons", outDir=opts.odir+"/truthAna")
    
    truth_is_match, reco_is_match = truth_reco_match(truth_electrons, signal_electrons)
    evt_mask = ak.num(truth_electrons)==2
    dMs = {"dm1" : (dM > 0.99) & (dM < 1.01), "dm3" : (dM > 2.99) & (dM < 3.01) }
    for dMname in dMs:
        dm_mask = dMs[dMname]
        pt1 = ak.max(truth_electrons[evt_mask & dm_mask].pt,axis=1)
        pt2 = ak.min(truth_electrons[evt_mask & dm_mask].pt,axis=1)
        # pt1_match = ak.max(truth_electrons[evt_mask & dm_mask].pt,axis=1)
        # pt2_match = ak.min(truth_electrons[evt_mask & dm_mask].pt,axis=1)
        plotHist("pt1", vals=pt1, xtitle="leading truth electron pt", outDir=opts.odir+"/truthAna/"+dMname)
        plotHist("pt2", vals=pt2, xtitle="subleading truth electron pt", outDir=opts.odir+"/truthAna/"+dMname)

        # matched_truth = truth_electrons[matching_truth_mask]
        # unmatched_truth = truth_electrons[~matching_truth_mask]
        # matched_reco = signal_electrons[matching_reco_mask]
        # unmatched_reco = signal_electrons[~matching_reco_mask]
        
    # STUDY RECO Working Points
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
        if opts.ROC: #doROC  cut_name in roc_config:
            maxeff = ak.count(truth_electrons)
            if maxeff: maxeff = ak.count(matched_truth) / maxeff
            rocs={}
            for var in roc_config:
                signal = ak.to_numpy(ak.flatten(matched_reco[var]))
                background = ak.to_numpy(ak.flatten(unmatched_reco[var]))
                roc = plotROC(var, signal, background, maxeff=maxeff, nEvts=nEvents, signalLike=roc_config[var],
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


        plotHist("met", vals=met, nbins=80, xtitle="met", outDir=opts.odir+"/met")
        plotHist("metPhi", vals=metPhi, nbins=80, xtitle="met phi", outDir=opts.odir+"/met")

        for dMname in dMs:
            dm_mask = dMs[dMname]
            # pt1 = ak.max(truth_electrons[evt_mask & dm_mask].pt,axis=1)
            
            nMatch = plotHist("n_matchedSignal_2eTruth_"+dMname+"_"+cut_name, vals=ak.num(matched_reco[evt_mask & dm_mask]),
                              lims=(-0.5,4.5), nbins=5,outDir=opts.odir+"/event", norm=True,
                              xtitle="matched signal electron multiplicity for "+dMname+" 2e events")
            print("Acceptance for 2e events for {} with '{}' selection is: {:.3g}".format(dMname,cut_name, nMatch[0][0][2]) )
        # plotHist("n_matchedSignal_"+cut_name, vals=ak.num(matched_reco), lims=(-0.5,4.5), nbins=5,
        #          xtitle="matched signal electron multiplicity", outDir=opts.odir+"/event")
        
    
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





    
    #combinePDFs()
    
    
    
    
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
    parser.add_option("--ROC", action='store_true', default = False, help="do ROCs only")
    parser.add_option("--drawMatchedTruth", action='store_true', default = False, help="histogram all matched/unmatched truth eles")
    (options, args) = parser.parse_args()
    analyze(options, args)
