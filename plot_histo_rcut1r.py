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


ifn0=['14','73']
ifn1=['g1','g3']
label0=['N=14','N=73']
label1=['(A)','(B)']
label2=[r'$L_{99}$',r'$L_{166}$']

nbin=50

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
            
def get_histo1(dist,density0):
    ''' density0= True: make PDF, 
        density0=False: histogram divided by nrun '''
    dmax=np.zeros((2,2))
    ndata=np.zeros((2,2),dtype=int)
    hist=[ [[] for k in range(2) ] for i in range(2)]
    bin_orig=[ [[] for k in range(2) ] for i in range(2)]
    edges=[ [[] for k in range(2) ] for i in range(2)]
    for i in range(len(dist)):
        for k in range(len(dist[i])):
            ndata[i][k]=len(dist[i][k])
            dmax[i][k]=np.amax(dist[i][k][:,1])
            bin_orig[i][k]=np.linspace(-0.5,dmax[i][k]+0.5,int(dmax[i][k])+2)
            hist[i][k],edges[i][k] = np.histogram(dist[i][k][:,1],bin_orig[i][k],density=density0)

    # if density0=False, hist: int64. Make a float array
    histf=[[] for k in range(2)]
    
    for i in range(len(dist)):
        histf[i]=[ [] for k in range(len(dist[i]))]
        for k in range(len(dist[i])):
            for p in range(len(hist[i][k])):
                if density0 == False:
                    rdum=float(hist[i][k][p])/float(len(dist[i][k]))
                else:
                    rdum=hist[i][k][p]
                histf[i][k].append(rdum)

    return ndata,dmax,histf,edges


def get_stat(hist,edges):
    avg=np.zeros((2,2)); std=np.zeros((2,2))
    for i in range(len(ifn0)):
        for k in range(len(ifn1)):
            avg0=0.; std0=0.; ptot=0.; # ptot: for normalization
            for p in range(len(hist[i][k])):
                p0=(edges[i][k][p+1]-edges[i][k][p])*hist[i][k][p] # probab
                ptot=ptot+p0 
                n0=0.5*(edges[i][k][p]+edges[i][k][p+1])
                avg0=avg0+n0*p0; std0=std0+n0*n0*p0
            avg[i][k]=avg0/ptot
            std0=std0/ptot-avg0*avg0; std[i][k]=np.sqrt(std0)
    return avg,std

def prep_fig(hist,edges):

    # set fig size proportional to number of rows & columns (WxH)
    fs0=20; fs1=16
    #bw0=0.3 # bar width

    fig = plt.figure(figsize=(8,4))
    gs = GridSpec(1,2)
    ax= []; ax2=[]
    ax.append(fig.add_subplot(gs[0,0]))
    ax.append(fig.add_subplot(gs[0,1]))

    #for i in range(len(hist[1][0])):
    #    print(edges[1][0][i],hist[1][0][i])
    #quit()
    
    for i in range(len(hist)):
        for k in range(2):
            ax[i].stairs(hist[i][k],edges[i][k],label=label2[k],linewidth=2)
        ax[i].set_xlabel(r'Number of pairs',fontsize=fs0)
        ax[i].tick_params(axis='both', which='major', labelsize=fs1)
    ax[0].legend(fontsize=fs0,frameon=False,loc=[0.5,0.08]) # ncol=2,1)
    ax[0].set_ylabel(r'Normalized Distribution',fontsize=fs0)

    for i in range(2):
        ax[i].text(0.72,0.91,label0[i],fontsize=fs0,transform=ax[i].transAxes)
        #ax[i].text(0.65,0.93,label1[i]+' '+label0[i],fontsize=fs0,transform=ax[i].transAxes)

        #alpha=0.15,ec='k',color='b',linewidth=2,capsize=6)

    fig.subplots_adjust(\
                        left  = 0.118 ,  # left margin
                        right = 0.995  ,  # start of right margin
                        bottom = 0.165 ,  # bottom margin
                        top = 0.992    ,  # start of top margin
                        wspace = 0.27 ,  # hgap
                        #hspace = 0.2   # vgap
    )
    #fig.tight_layout() # remove margins
    return fig

def read_data():
    dist=[[] for i in range(2)]
    for i in range(len(ifn0)):
        sdum0=ifn0[i]
        for sdum1 in ifn1:
            sdum='./data1/out'+sdum0+sdum1+'_rcut.dat'
            a=np.loadtxt(sdum)
            dist[i].append(a)
    return dist

def write_histo(ndata,dmax,hist,edges,avg,std):
    for i in range(len(ifn0)):
        sdum0=ifn0[i]
        for k in range(len(ifn1)):
            sdum1=ifn1[k]
            sdum='./data1/histo'+sdum0+sdum1+'_rcut.dat'
            ff = open(sdum,'w')
            ff.write('# Number of data points: {:d}\n'.format(ndata[i][k]))
            ff.write('# maxval= {:.5f}  nbin= {:d}\n'.format(dmax[i][k],nbin))
            ff.write('# rcut= '+label2[i]+'\n')
            ff.write('# Edges: min= {:.8f}  max= {:.8f}\n'.format(edges[i][k][0],edges[i][k][-1]))
            ff.write('# avg/std of counts: {:.5f} {:.5f}\n'.format(avg[i][k],std[i][k]))
            ff.write('#count  histo(normalized)\n')
            norm0=0.;
            for p in range(len(hist[i][k])):
                norm0=norm0+(edges[i][k][p+1]-edges[i][k][p])*hist[i][k][p]
                # len(edges[i]) is 1 larger than len(hist[i])
                rdum=0.5*(edges[i][k][p]+edges[i][k][p+1]) 
                sdum1='{0:5.1f}  {1:.5e}\n'.format(rdum,hist[i][k][p])
                ff.write(sdum1)
            sdum2=': cumulative probab= {:f}'.format(norm0)
            print(sdum0+ifn1[k]+sdum2)
    
#######################################################
if __name__=="__main__":
    #idir, pull_ini,pull_fin = get_args()
    dist = read_data()
    ndata,dmax,hist,edges = get_histo1(dist,False)
    avg,std = get_stat(hist,edges)
    write_histo(ndata,dmax,hist,edges,avg,std)

    fig1=prep_fig(hist,edges)
    fig1.savefig('histo_rcut1r.pdf',format='pdf')
    quit()    
