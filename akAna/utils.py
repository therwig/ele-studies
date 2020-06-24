import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import awkward1 as ak
from scipy.stats import beta
from math import ceil

def plot(vals,
         savename,
         nbins=40,
         xtitle='',
         ytitle='Entries',
         outDir='plots',
         formats=['pdf']
         ):
    plt.figure(figsize=(6,4))
    plt.hist(vals, nbins)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    for form in formats:
        sname="{}/{}.{}".format(outDir,savename,form)
        plt.savefig(sname)
        print("Saved: "+sname)
    plt.close()

def plotCollection(objs,
                   savename,
                   xtitle='',
                   ):
    plot(ak.num(objs),"n_"+savename, xtitle=xtitle+" multiplicity")
    for attr in objs.columns:
        plot(ak.flatten(objs[attr]), savename+'_'+attr, xtitle=xtitle+" "+attr)
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

def efficiency(passVals, totVals, savename,
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
    for form in formats:
        sname="{}/{}.{}".format(outDir,savename,form)
    plt.savefig(sname)
    print("Saved: "+sname)
    plt.close()
B
