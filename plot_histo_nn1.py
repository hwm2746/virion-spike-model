''' 
 This file is part of virion-spike-model, which is distributed under the
 MIT license, as described in the LICENSE file in the top level directory
 of this project.
 
 Author: Wonmuk Hwang

'''
import sys
import argparse
import numpy as np
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.gridspec import GridSpec


ifn=['14','73'] #ifn=['14a','73a']
label0=['N=14','N=73']
label1=['(A)','(B)']

nbin=50

#TOL=1.0e-8 # tolerance for setting probability equal to zero
#dbin=[1,2,4,6,8,10,12,15]
#ndihe=6 # number of dihedral angles

#def get_args():
#    """
#    Process the command line options.
#    """
#
#    parser = argparse.ArgumentParser(description=
#               'Usage: histogram.py -i ifname -o ofname -b dbin')
#    parser.add_argument('-i', metavar='idir')
#    parser.add_argument('-pi', metavar='pull_ini')
#    parser.add_argument('-pf', metavar='pull_fin')
#
#    args = parser.parse_args() # args is a dictionary
#    #    dbin = float(vars(args)['b'])
#    idir = vars(args)['i']
#    pull_ini = int(vars(args)['pi'])
#    pull_fin = int(vars(args)['pf'])
#
#    return idir, pull_ini, pull_fin

def get_histo(dist,nbin):
    dmax=np.zeros(2)
    ndata=np.zeros(2,dtype=int)
    hist=[[] for i in range(2)]
    bin_orig=[ [] for i in range(2)]
    edges=[[] for i in range(2)]
    for i in range(len(dist)):
        ndata[i]=len(dist[i][:,1])
        dmax[i]=np.amax(dist[i][:,1])
        rdum=0.5*dmax[i]/float(nbin) # pad for min/max
        bin_orig[i]=np.linspace(0.-rdum,dmax[i]+rdum,nbin)
        #print(bin_orig[i][0],bin_orig[i][-1])
        hist[i],edges[i] = np.histogram(dist[i][:,1],bin_orig[i],density=True)
        #print(edges[i])
        #print(np.amax(dist[i][:,1]))

    hcumul=np.zeros_like(hist) # cumulative sum       
    for i in range(len(hist)): 
        rdum=0.;
        for k in range(len(hist[i])):
            rdum=rdum+hist[i][k]*(edges[i][k+1]-edges[i][k])
            hcumul[i][k]=rdum

    return ndata,dmax,hist,edges,hcumul

def prep_fig(hist,edges,hcumul):

    # set fig size proportional to number of rows & columns (WxH)
    fs0=14; fs1=12
    #bw0=0.3 # bar width

    fig = plt.figure(figsize=(8,4))
    gs = GridSpec(1,2)
    ax= []; ax2=[]
    ax.append(fig.add_subplot(gs[0,0]))
    ax.append(fig.add_subplot(gs[0,1]))

    e0=[ [] for i in range(len(hist))]
    for i in range(len(hist)):
        #e0[i]=np.zeros(len(hist[i]))
        #for k in range(len(e0)):
        #    e0[i][k]=0.5*(edges[i][k]+edges[i][k+1]) 
        ax[i].stairs(hist[i],edges[i],color='b')
        #ax[i].legend(fontsize=fs1,frameon=False) #,loc=[0.27,0.6]) # ncol=2,1)

        ax[i].set_xlabel(r'NN Dist (nm)',fontsize=fs0)
        #alpha=0.15,ec='k',color='b',linewidth=2,capsize=6)

    for i in range(len(hist)): # plot cumulative sum
        # twin object for two different y-axis on the sample plot
        ax2.append(ax[i].twinx())
        # make a plot with different y-axis using second axis object
        ax2[i].stairs(hcumul[i],edges[i],color="r")

    for i in range(len(hist)):
        ax[i].tick_params(axis='y', colors='b')
        ax2[i].tick_params(axis='y', colors='r')
        ax[i].text(0.7,0.8,label0[i],fontsize=fs0,transform=ax[i].transAxes)
        ax[i].text(-0.1,-.12,label1[i],fontsize=fs0,transform=ax[i].transAxes)

    ax[0].set_ylabel(r'Normalized Distribution',color='b',fontsize=fs0)
    ax2[1].set_ylabel('Cumul. Distribution',color='r',fontsize=fs0)
    

#
#    # custom x-axis label
#    ax1.xaxis.set_tick_params(labelsize=fs0) # set xtick label size
#    plt.xticks(x0,mname0) # fontsize doesn't work here
#    
#    ax1.yaxis.set_tick_params(labelsize=fs0)
#    ax2.yaxis.set_tick_params(labelsize=fs0)
#    #ax2.set_ylim(0.,3.)
#    #
#    #
#    ax1.set_ylabel(r'V$\beta$-FG Contact Count',fontsize=fs0,color='b')
#    ax2.set_ylabel(r'C$\beta$-FG Contact Count',fontsize=fs0,color='m')
#
    fig.subplots_adjust(\
                        left  = 0.1 ,  # left margin
                        right = 0.93  ,  # start of right margin
                        bottom = 0.12 ,  # bottom margin
                        top = 0.985    ,  # start of top margin
                        wspace = 0.3 ,  # hgap
                        #hspace = 0.2   # vgap
    )
    #fig.tight_layout() # remove margins
    return fig

def read_data():
    dist=[]
    for ifn0 in ifn:
        sdum='./data1/out'+ifn0+'g1_dist.dat'
        a=np.loadtxt(sdum)
        dist.append(a)
    return dist

def write_histo(ndata,dmax,hist,edges,hcumul):
    for i in range(len(ifn)):
        ifn0=ifn[i]
        sdum='./data1/histo'+ifn0+'.dat'
        ff = open(sdum,'w')
        ff.write('# Number of data points: {:d}\n'.format(ndata[i]))
        ff.write('# maxval= {:.5f}  nbin= {:d}\n'.format(dmax[i],nbin))
        ff.write('# Edges: min= {:.8f}  max= {:.8f}\n'.format(edges[i][0],\
                                                              edges[i][-1]))
        ff.write('#dist(nm)  histo(normalized) cumulative_histo (dist: bin center)\n')
        for k in range(len(hist[i])):
            # len(edges[i]) is 1 larger than len(hist[i])
            rdum=0.5*(edges[i][k]+edges[i][k+1]) 
            sdum1='{0:11.8f}  {1:.5e} {2:.5e}\n'.format(rdum,hist[i][k],
                                                      hcumul[i][k])
            ff.write(sdum1)
    
#######################################################
if __name__=="__main__":
    #idir, pull_ini,pull_fin = get_args()
    dist = read_data()
    ndata,dmax,hist,edges,hcumul = get_histo(dist,nbin)
    write_histo(ndata,dmax,hist,edges,hcumul)
    fig1=prep_fig(hist,edges,hcumul)
    #fig1.savefig('histo_nn.png',format='png')
    fig1.savefig('histo_nn1.pdf',format='pdf')
    quit()
    
    #print np.sum(hist) # shows that hist is PDF, not discrete probability
    #write_histo(hist,bin_edges,ofname)

    #for i in range(len(dbin)):
    #    hist, bin_edges = get_histo(ifname,dbin[i])
    #    get_entropy(hist,dbin[i])
    #entropy0 = estimate_entropy(ifname,dbin)
