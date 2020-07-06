import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import awkward1 as ak
import scipy
from scipy.stats import beta
from math import ceil
import pathlib
import os

from plot_config import config,notNone

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
        if notNone(config[var].nbins): nbins = config[var].nbins
        if notNone(config[var].lo) and notNone(config[var].hi): lims = (config[var].lo, config[var].hi)
        if notNone(config[var].log): isLog = config[var].log
        if notNone(config[var].name): xtitle = config[var].name
        #print('for var '+var+' setting',config[var].nbins,config[var].lo, config[var].hi)

    # plot and postfix overrides
    if vals:
        packed = plt.hist(vals, nbins, range=lims, log=isLog, density=norm, histtype='step', label=leg)
        print('doing',nbins, lims)
    elif hists:
        for packed_hist in hists:
            vals, bins, patches = packed_hist
            centers = (bins[1:] + bins[:-1])/2
            plt.hist(x=centers, weights=vals, bins=bins, log=isLog, density=norm, histtype='step', label=leg)
            # plt.hist(x=np.ones_like(vals), weights=vals, bins=bins, log=isLog)
            #plt.hist(data=vals, bins=bins, log=isLog)
        packed=None
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)

    sname="{}/{}".format(outDir,savename)
    sname = sname.replace('//','/')
    if writeTitle: plt.text(-0.1, 1.1, sname, fontsize=6, transform=plt.gca().transAxes)
    
    if leg: fig.legend()
    
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

def plotProfile(savename,
                xs=None,
                ys=None,
                hists=None,
                nbins=15,
                var='',
                collection='',
                xtitle='',
                lims=(0,30),
                isLog=False,
                ytitle='Entries',
                outDir='plots',
                leg=None,
                writeTitle=True,
                formats=['pdf']
):

    fig = plt.figure(figsize=(6,4))
    if var in config:
        if config[var].nbins: nbins = config[var].nbins
        if config[var].lo and config[var].hi: lims = (config[var].lo, config[var].hi)
        if config[var].log: isLog = config[var].log
        if config[var].name: ytitle = config[var].name

    profs=[]
    for i,x in enumerate(xs):
        #if i>0: continue
        y=ys[i]
        if lims==None: lims = (np.min(x),np.max(x))
        # if lims==None: lims = (x.min(),x.max())
        median_result = scipy.stats.binned_statistic(x, y, bins=nbins, range=lims, statistic=lambda x: np.quantile(x,0.5))
        lo_result     = scipy.stats.binned_statistic(x, y, bins=nbins, range=lims, statistic=lambda x: np.quantile(x,0.5-0.68/2))
        hi_result     = scipy.stats.binned_statistic(x, y, bins=nbins, range=lims, statistic=lambda x: np.quantile(x,0.5+0.68/2))
        median = np.nan_to_num(median_result.statistic, posinf=0, neginf=0)
        hi = np.nan_to_num(hi_result.statistic, posinf=0, neginf=0)
        lo = np.nan_to_num(lo_result.statistic, posinf=0, neginf=0)
        hie = hi-median
        loe = median-lo
        bin_edges = median_result.bin_edges
        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.

        offset = 0.33*(bin_centers[1]-bin_centers[0])
        off = offset * (i-1)/2 * (-1. if i%2 else 1.)
        plt.errorbar(x=bin_centers+off, y=median, yerr=[loe,hie], label=leg[i]) #, linestyle='none', marker='.', label=leg)
        
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)    
    if leg: fig.legend()

    sname="{}/{}".format(outDir,savename)
    sname = sname.replace('//','/')
    if writeTitle: plt.text(-0.1, 1.1, sname, fontsize=6, transform=plt.gca().transAxes)

    # save figs
    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
    global pdflist
    for form in formats:
        snames=sname+"."+form
        plt.savefig(snames)
        print("Saved: "+snames)
        if form=='pdf': pdflist.append(snames)
    plt.close()
    return


def plotCollection(objs,
                   savename,
                   xtitle='',
                   leg=None,
                   outDir='plots',
                   normAttrs=False,
                   plotlist=None,
                   profile=False,
                  ):
    if not type(objs) is list:
        objs = [objs]
        
    plotHist("n_"+savename, vals=[ak.num(o) for o in objs], xtitle=xtitle+" multiplicity", leg=leg, outDir=outDir)
    if plotlist is None:
        plotlist = {x:True for x in objs[0].columns}
    for attr in plotlist:
        pltFunc = plotlist[attr]
        if type(pltFunc) is bool:
            if not pltFunc: continue
            pltFunc = lambda x : x[attr]
        plotHist(savename+'_'+attr, var=attr, vals=[ak.to_list(ak.flatten(pltFunc(o))) for o in objs], xtitle=xtitle+" "+attr, leg=leg, outDir=outDir, norm=normAttrs)
        #break # to just plot one thing
        if profile and (attr in config) and (config[attr].profile):
            plotProfile(savename+'_'+attr+'_prof', var=attr,
                        xs=[ak.flatten(o['pt']) for o in objs],
                        ys=[ak.flatten(pltFunc(o)) for o in objs],
                        xtitle=xtitle+" p_T [GeV]", leg=leg, outDir=outDir+"/profiles/")


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
        os.system(cmd)
        print("Combined plots from this latest run into " + outname + ".pdf")
