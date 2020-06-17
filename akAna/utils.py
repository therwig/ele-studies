import numpy as np
import matplotlib.pyplot as plt

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
        plt.savefig("{}/{}.{}".format(outDir,savename,form))
    plt.close()
