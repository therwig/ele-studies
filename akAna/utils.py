import numpy as np
import matplotlib.pyplot as plt
import awkward1 as ak

def plot(vals,
         savename,
         nbins=40,
         xtitle='',
         ytitle='Entries',
         outDir='plots',
         formats=['png']
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
