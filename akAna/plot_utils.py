import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import awkward1 as ak
from scipy.stats import beta
from math import ceil
import pathlib
import os

from plot_config import config

pdflist=[]

def plotHist(savename,
             vals=None,
             hists=None,
             nbins=40,
             var='',
             collection='',
             xtitle='',
             lims=None,
             isLog=False,
             norm=False,
             ytitle='Entries',
             outDir='plots',
             leg=None,
             writeTitle=True,
             formats=['pdf']
            ):

    fig = plt.figure(figsize=(6,4))

    if var in config:
        nbins = config[var].nbins
        lims = (config[var].lo, config[var].hi)
        isLog = config[var].log
        xtitle = config[var].name

    # plot and postfix overrides
    if vals:
        packed = plt.hist(vals, nbins, range=lims, log=isLog, density=norm)
    elif hists:
        for packed_hist in hists:
            vals, bins, patches = packed_hist
            centers = (bins[1:] + bins[:-1])/2
            plt.hist(x=centers, weights=vals, bins=bins, log=isLog, density=norm)
            # plt.hist(x=np.ones_like(vals), weights=vals, bins=bins, log=isLog)
            #plt.hist(data=vals, bins=bins, log=isLog)
        packed=None
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)

    sname="{}/{}".format(outDir,savename)
    sname = sname.replace('//','/')
    if writeTitle: plt.text(-0.1, 1.1, sname, fontsize=6, transform=plt.gca().transAxes)
    
    if leg: fig.legend(leg)
    
    # save figs
    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
    global pdflist
    for form in formats:
        snames=sname+"."+form
        plt.savefig(snames)
        print("Saved: "+snames)
        if form=='pdf': pdflist.append(snames)
    plt.close()
    return packed

def plotCollection(objs,
                   savename,
                   xtitle='',
                   leg=None,
                   outDir='plots',
                   normAttrs=False,
                  ):
    if type(objs) is list:
        plotHist("n_"+savename, vals=[ak.num(o) for o in objs], xtitle=xtitle+" multiplicity", leg=leg, outDir=outDir)
        for attr in objs[0].columns:
            plotHist(savename+'_'+attr, vals=[ak.to_list(ak.flatten(o[attr])) for o in objs], xtitle=xtitle+" "+attr, leg=leg,
                     outDir=outDir, norm=normAttrs)
    else:
        plotHist("n_"+savename, vals=ak.num(objs), xtitle=xtitle+" multiplicity", leg=leg, outDir=outDir)
        for attr in objs.columns:
            plotHist(savename+'_'+attr, vals=ak.flatten(objs[attr]), xtitle=xtitle+" "+attr, leg=leg,
                     outDir=outDir, norm=normAttrs)
        #break # to just plot one thing

def ErrDivide(p,n, lvl=0.68):
    if n < p: raise Exception("ErrDivide: Pass greater than total!")
    if n<0 or p<0: raise Exception("ErrDivide: Negative yields!")
    r = p/n
    a=p+1 # p + alpha
    b=(n-p)+1 # n-p + bets
    eLo = beta.ppf((1-lvl)/2, a, b)
    eHi = beta.ppf((1+lvl)/2, a, b)
    eLoRnd = int(eLo*n)/n
    eHiRnd = ceil(eHi*n)/n
    return r, eLoRnd, eHiRnd #eLo, eHi

def plotEfficiency(savename,
                   effs=None,
                   passVals=None, totVals=None,
                   nbins=10, lims=None,
                   xtitle='',
                   outDir='plots',
                   leg=None,
                   formats=['pdf']):

    fig = plt.figure(figsize=(6,4))
    
    if passVals and totVals:
        if lims==None: lims = (totVals.min(),totVals.max())
        bin_edges = np.linspace(lims[0],lims[1],nbins+1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
        passHist = np.histogram(passVals,bin_edges)
        totHist = np.histogram(totVals,bin_edges)
        effs = np.array([ErrDivide(passHist[0][i],totHist[0][i]) for i in range(len(bin_centers))])
        eff, eDn, eUp = np.hsplit(effs,3)
        plt.errorbar(bin_centers, eff[:,0], yerr=[(eff-eDn)[:,0], (eUp-eff)[:,0]])
        packed = [(bin_centers, (eff[:,0], (eff-eDn)[:,0], (eUp-eff)[:,0]) )]
    elif effs:
        for packed_eff in effs:
            bins, eff = packed_eff
            plt.errorbar(bins, eff[0], yerr=[eff[1],eff[2]])
        packed = effs
    else:
        raise Exception("Must pass either pass/total vals or existing effs!")
    
    if leg: fig.legend(leg)
    plt.xlabel(xtitle)
    plt.ylabel('Efficiency')
    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
    global pdflist
    for form in formats:
        sname="{}/{}.{}".format(outDir,savename,form)
        sname = sname.replace('//','/')
        plt.savefig(sname)
        print("Saved: "+sname)
        if form=='pdf': pdflist.append(sname)
    plt.close()

    return packed

def combinePDFs(outname='latest'):
    '''
    Run like:
    pdfjoin output_plots/*pdf -o out2.pdf
    pdfunite *pdf outUNITE.pdf
    '''
    global pdflist
    inlist = " ".join(pdflist)
    user = os.environ['USER']
    #host = os.environ['HOSTNAME']
    if 'herwig' in user:
        cmd = "pdfjoin -q "+inlist+" -o "+outname + ".pdf"
        if len(pdflist) >= 248:
            print("Can't combine more than 248 inputs!")
        else:
            os.system(cmd)
            print("Combined plots from this latest run into " + outname + ".pdf")
    else:
        cmd = "pdfunite "+inlist+" "+outname + ".pdf"
        print("Combined plots from this latest run into " + outname + ".pdf")
