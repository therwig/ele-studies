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
            n_objs = ak.num(cms_events[collection]['pt'])
            plotHist(n_objs,"n_"+collection, xtitle=collection+" multiplicity", outDir=opts.odir)
            for attr in cms_events[collection].columns:
                flat_values = ak.flatten(ak.to_list(cms_events[collection][attr]))
                plotHist(flat_values,collection+'_'+attr, var=attr, xtitle=collection+" "+attr, outDir=opts.odir)
    
    
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
        plotHist(ak.num(truth_electrons),"n_genElectrons", xtitle="gen electron multiplicity", outDir=opts.odir)
        for attr in truth_electrons.columns:
            plotHist(ak.flatten(ak.to_list(truth_electrons[attr])),'genElectron_'+attr, var=attr, xtitle="Gen Electron "+attr, outDir=opts.odir)
    
    # DEFINE THE RECO ELECTRONS        
    #reco_electrons = cms_events['electrons']
    reco_electrons = cms_events['softElectrons']
    
    if opts.drawRecoElectrons:
        plotHist(ak.num(reco_electrons), "n_recoElectrons", xtitle="reco electron multiplicity", outDir=opts.odir)
        for recoe in reco_electrons.columns:
            plotHist(ak.flatten(ak.to_list(reco_electrons[recoe])),'recoElectron_'+recoe, var=recoe, xtitle="Reco Electron"+recoe, outDir=opts.odir)
    
    # BEGIN THE ANALYSIS
    evt_mask = ak.num(truth_electrons)==2
    
    match_builder = dr_match(truth_electrons, reco_electrons, ak.ArrayBuilder())
    match_extension = match_builder.snapshot()
    matching_truth_mask = match_extension.reco_index >= 0
    matching_reco_mask = match_extension.reco_index[matching_truth_mask]
    
    matched_truth = truth_electrons[matching_truth_mask]
    unmatched_truth = truth_electrons[~matching_truth_mask]
    matched_reco = reco_electrons[matching_reco_mask]
    
    if opts.drawMatchedCollections:
        plotCollection(matched_truth,   "matched_truth_ele",   xtitle="matched truth electron")
        plotCollection(unmatched_truth, "unmatched_truth_ele", xtitle="unmatched truth electron")
        plotCollection(matched_reco,    "matched_truth_reco",  xtitle="matched reco electron")
    
    # display matching efficiencies
    passVals = ak.to_list(ak.flatten(matched_truth.pt))
    totVals = ak.to_list(ak.flatten(truth_electrons.pt))
    plotEfficiency(passVals, totVals, "eff_pt", lims=(0,20), nbins=20, xtitle="truth electron p_T [GeV]", outDir=opts.odir)

    combinePDFs()



    
if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-i',"--input", type="string", default = '', help="path to input file")
    parser.add_option('-o',"--odir", type="string", default = 'output_plots/', help="plots directory")
    parser.add_option("--drawInputHistograms", action='store_true', default = False, help="histogram ALL input collections")
    parser.add_option("--drawTruthElectrons", action='store_true', default = False, help="histogram selected truth electrons")
    parser.add_option("--drawRecoElectrons", action='store_true', default = False, help="histogram selected reco electrons")
    parser.add_option("--drawMatchedCollections", action='store_true', default = False, help="histogram all matched/unmatched truth and reco eles")
    (options, args) = parser.parse_args()
    analyze(options, args)
    
