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

def plotHist(vals,
             savename,
             nbins=40,
             var='',
             collection='',
             xtitle='',
             lims=None,
             isLog=False,
             ytitle='Entries',
             outDir='plots',
             formats=['pdf']
            ):

    plt.figure(figsize=(6,4))

    if var in config:
        nbins = config[var].nbins
        lims = (config[var].lo, config[var].hi)
        isLog = config[var].log
        xtitle = config[var].name

    # plot and postfix overrides
    plt.hist(vals, nbins, range=lims, log=isLog)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)

    # save figs
    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
    global pdflist
    for form in formats:
        sname="{}/{}.{}".format(outDir,savename,form)
        sname = sname.replace('//','/')
        plt.savefig(sname)
        print("Saved: "+sname)
        if form=='pdf': pdflist.append(sname)
    plt.close()

def plotCollection(objs,
                   savename,
                   xtitle='',
                   ):
    plot(ak.num(objs),"n_"+savename, xtitle=xtitle+" multiplicity")
    for attr in objs.columns:
        plotHist(ak.flatten(objs[attr]), savename+'_'+attr, xtitle=xtitle+" "+attr)
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

def plotEfficiency(passVals, totVals, savename,
               nbins=10, lims=None,
               xtitle='',
               outDir='plots',
               formats=['pdf']):
    if lims==None: lims = (totVals.min(),totVals.max())
    bin_edges = np.linspace(lims[0],lims[1],nbins+1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
    passHist = np.histogram(passVals,bin_edges)
    totHist = np.histogram(totVals,bin_edges)
    effs = np.array([ErrDivide(passHist[0][i],totHist[0][i]) for i in range(len(bin_centers))])
    eff, eDn, eUp = np.hsplit(effs,3)
    
    plt.figure(figsize=(6,4))
    plt.errorbar(bin_centers, eff[:,0], yerr=[(eff-eDn)[:,0], (eUp-eff)[:,0]])
    
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
    else:
        cmd = "pdfunite "+inlist+" "+outname + ".pdf"
    os.system(cmd)
    print("Combined plots from this latest run into " + outname + ".pdf")
