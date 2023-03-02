#!/usr/bin/python2.7


import sys,math,os,string
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import random
import numpy as np


# filename = sys.argv[1]

fig,ax=plt.subplots(3,1)

for filecounter,filename in enumerate(sys.argv[1:]):
    time = 0
    la=[]
    l_alive=[]
    lr=[]
    lR=[]
    ldistrR=[]
    ltime=[]
    with open(filename,"r") as fin:
        print "Opening file: ", filename
        # this reads one line at a time
        for line in fin:
            line = line.split()
            timenow=int(line[0])
            
            if timenow != time:
                time = timenow
                ltime.append(time)
                l_alive.append( np.mean(la) )
                lR.append( np.mean(lr) )
                ldistrR.append(lr)
                lr=[]
                la=[]
            la.append( 1 if float(line[1])>0.0000001 else 0 )
            lr.append( float(line[1]) )
            
    # print l_exp_to_mut
    linestyle='-'
    # print ltime
    # print l_alive
    ax[0].plot(ltime, l_alive)
    ax[0].plot(ltime, lR,label=filename)
    ax[0].set_ylim([0,1])
    ax[0].boxplot(ldistrR,positions=ltime, showfliers=True,widths=(ltime[1]-ltime[0])/3.)
    # ax[0].plot(ltime, l_fitness,label=anc_string+"av fitness "+filename,color='royalblue',linestyle=linestyle)
    # if not "ancestor_trace" in filename:
    # 
    #     ax[0].plot(ltime, [a+b for a,b in zip(l_frac_expr,l_std_frexpr)],color='firebrick')
    #     ax[0].plot(ltime, [a-b for a,b in zip(l_frac_expr,l_std_frexpr)],color='firebrick')
    # 
    #     ax[0].plot(ltime, [a+b for a,b in zip(l_fitness,l_std_fitness)],color='royalblue')
    #     ax[0].plot(ltime, [a-b for a,b in zip(l_fitness,l_std_fitness)],color='royalblue')
    # 
    # ax[0].plot(ltime, l_open, label=anc_string+"func, open fract",color='darkgoldenrod',linestyle=linestyle)
    # ax[0].plot(ltime, l_bias, label=anc_string+"bias", color='seagreen',linestyle=linestyle)
    # ax[0].plot(ltime, l_exp_to_mut, label=anc_string+"exp_to_mut", color='magenta',linestyle=linestyle)
    # 
    ax[0].legend()
    # 
    # color = (random.random(),random.random(),random.random())
    # 
    # ax[1].plot(ltime,lav_a,label=anc_string+"av # func "+filename,linestyle='--',c=color)
    # ax[1].plot(ltime,lav_i,label=anc_string+"av # brok "+filename,linestyle='-',c=color )
    # if( len(sys.argv[1:]) < 4 and not  "ancestor_trace" in filename ):
    #     ax[1].plot(ltime,[a+b for a,b in zip(lav_a,l_std_a)],linestyle='--',c=color)
    #     ax[1].plot(ltime,[a-b for a,b in zip(lav_a,l_std_a)],linestyle='--',c=color)
    #     ax[1].plot(ltime,[a+b for a,b in zip(lav_i,l_std_i)],linestyle='-',c=color )
    #     ax[1].plot(ltime,[a-b for a,b in zip(lav_i,l_std_i)],linestyle='-',c=color )
    # 
    # ax[1].legend()
    # if 'ancestor_trace' in filename:
    #     # print lgenetag
    #     # timestep=ltime[1]-ltime[0]
    #     # print lgenetag
    #     max_etm = max([max(x) for x in lgenetag])
    #     min_etm = min([min(x) for x in lgenetag])
    #     print "min_etm,max_etm",min_etm,max_etm
    #     hist_2d = []
    #     for lgenetag_x in lgenetag:
    #         hist,bins = np.histogram(lgenetag_x, bins=np.linspace(min_etm,max_etm,100.))
    #         hist_2d.append(hist)
    #     # transpose 
    #     hist_2d = zip(*hist_2d)
    #     hist_2d = np.ma.masked_equal(hist_2d,0)
    #     print "hist bins[:-1] length = " , len(bins[:-1])
    #     print "ltime length" ,len(ltime)
    #     print "len(hist_2d)",len(hist_2d)
    #     print "len(hist_2d[0])",len(hist_2d[0])
    #     x, y = np.meshgrid(ltime,bins[:-1])
    #     ax[2].pcolormesh(x,y,hist_2d)#,label='anc. neutr. mark')
    #     # ax[2].text((max_etm+min_etm)/2 - (max_etm-min_etm)/2, 1, 'anc. neutr. mark', bbox={'facecolor': 'white', 'pad': 10})
    ax[0].margins(0,0)
    ax[1].margins(0,0)
    ax[2].margins(0,0)
plt.show()
