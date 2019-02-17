#!/usr/bin/env python

import os
import sys
import sys, os, time
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from astropy import units as u
from generic_input_pars import Generic_Pars
from build_fsps_model import Build_Fsps

sys.path.append(
  os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from Dcolor2sSNRL_gen import Generate_Curve
from SN_rate import Model_Rates

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'  
fs = 24.
lw = 3.   
taus = [1., 1.5, 2., 3., 4., 5., 7., 10.]

s1s2 = zip([-1., -1.5, -2., -1., -1.],[-1., -1.5, -.5, -2.,-3.])
label = [r'-1/-1', r'-1.5/-1.5', r'-2/-0.5', r'-1/-2', r'-1/-3']

cZ = ['#fc9272','#ef3b2c','#a50f15']
ct = ['#a6bddb','#3690c0','#016450']
lsZ = ['--', '-', '-.']
lst = ['--', '-.', '-']

class Plot_Tests(object):
    """
    Description:
    ------------
    TBW.

    Parameters:
    -----------
    show_fig : ~bool
        True of False. Whether to show or not the produced figure.
    save_fig : ~bool
        True of False. Whether to save or not the produced figure.

    Notes:
    ------
    The normalization used here is different than in Heringer+ 2017. In that
    paper the DTD is normalized at 0.5 Gyr, whereas here an arbitraty
    constant (A=10**-12) is given to the DTD of the form SN rate = A*t**s1
    where s1 is the slope prior to 1 Gyr. 
             
    Outputs:
    --------
    ./../OUTPUT_FILES/ANALYSES_FIGURES/Fig_sSNRL.pdf
    
    References:
    -----------
    Heringer+ 2017: http://adsabs.harvard.edu/abs/2017ApJ...834...15H
    """         
    def __init__(self, show_fig=True, save_fig=False):
        self.show_fig = show_fig
        self.save_fig = save_fig 
       
        self.fig, (self.ax1, self.ax2) = plt.subplots(
          1,2, figsize=(16,8), sharey=True)
      
        self.make_plot()
                
    def set_fig_frame(self):

        plt.subplots_adjust(wspace=0.03)
        
        x_label = r'$\Delta (g-r)$'
        y_label = r'$\rm{log\ sSNR}_L\ \rm{[yr^{-1}\ L_\odot^{-1}]}$'
        
        self.ax1.set_xlabel(x_label, fontsize=fs)
        self.ax1.set_ylabel(y_label, fontsize=fs)
        self.ax1.set_xlim(-1.05,0.35)
        self.ax1.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax1.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax1.tick_params('both', length=12, width=2., which='major',
                                 direction='in', right=True, top=True)
        self.ax1.tick_params('both', length=6, width=2., which='minor',
                                 direction='in', right=True, top=True) 
        self.ax1.xaxis.set_minor_locator(MultipleLocator(.05))
        self.ax1.xaxis.set_major_locator(MultipleLocator(.2))
        self.ax1.yaxis.set_minor_locator(MultipleLocator(.25))
        self.ax1.yaxis.set_major_locator(MultipleLocator(.5))  

        self.ax2.set_xlabel(x_label, fontsize=fs)
        self.ax2.set_xlim(-1.05,0.35)
        self.ax2.set_ylim(-14.75,-10.75)
        self.ax2.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax2.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax2.tick_params('both', length=12, width=2., which='major',
                                 direction='in', right=True, top=True)
        self.ax2.tick_params('both', length=6, width=2., which='minor',
                                 direction='in', right=True, top=True) 
        self.ax2.xaxis.set_minor_locator(MultipleLocator(.05))
        self.ax2.xaxis.set_major_locator(MultipleLocator(.2))
        self.ax2.yaxis.set_minor_locator(MultipleLocator(.25))
        self.ax2.yaxis.set_major_locator(MultipleLocator(.5))  
        plt.setp(self.ax2.get_yticklabels(), visible=False)

    def plot_models(self):

        for k, Z in enumerate(['0.0150', '0.0190', '0.0300']):
            _inputs = Generic_Pars('exponential', 'Kroupa', Z, 100.e6 * u.yr)
            _D = Build_Fsps(_inputs).D        

            for i, (s1,s2) in enumerate(s1s2):
                Sgen = Generate_Curve(_inputs, _D, s1, s2)
                x = Sgen.Dcd_fine
                y = np.log10(Sgen.sSNRL_fine * 1.e-12)
                self.ax1.plot(x, y, ls=lsZ[k], marker='None',
                              color=cZ[k], linewidth=lw, zorder=3)   
                if k==2:
                    self.ax1.text(0.05, y[-1] + 0.05, label[i], color='k', fontsize=fs)
        
        for k, t_ons in enumerate([40.e6 * u.yr, 70.e6 * u.yr, 100.e6 * u.yr]):
            _inputs = Generic_Pars('exponential', 'Kroupa', '0.0190', t_ons)
            _D = Build_Fsps(_inputs).D        

            for i, (s1,s2) in enumerate(s1s2):
                Sgen = Generate_Curve(_inputs, _D, s1, s2)
                x = Sgen.Dcd_fine
                y = np.log10(Sgen.sSNRL_fine * 1.e-12)
                self.ax2.plot(x, y, ls=lst[k], marker='None',
                              color=ct[k], linewidth=lw, zorder=3)   
                if k==2:
                    self.ax2.text(0.05, y[-1] + 0.05, label[i], color='k', fontsize=fs)

    def make_legend(self):
        self.ax1.plot(
          [np.nan], [np.nan], color=cZ[0], ls=lsZ[0], lw=lw,
          marker='None', label=r'$Z=0.015$')
        self.ax1.plot(
          [np.nan], [np.nan], color=cZ[1], ls=lsZ[1], lw=lw,
          marker='None', label=r'$Z=0.019$')
        self.ax1.plot(
          [np.nan], [np.nan], color=cZ[2], ls=lsZ[2], lw=lw,
          marker='None', label=r'$Z=0.030$')          
        self.ax1.legend(
          frameon=False, fontsize=fs, numpoints=1, ncol=1, labelspacing=.2,
          handlelength=1.5, handletextpad=.8, loc=3, bbox_to_anchor=(.2, 0.)) 

        self.ax2.plot(
          [np.nan], [np.nan], color=ct[0], ls=lst[0], lw=lw,
          marker='None', label=r'$t_{\rm{WD}}=40\, \rm{Myr}$')         

        self.ax2.plot(
          [np.nan], [np.nan], color=ct[1], ls=lst[1], lw=lw,
          marker='None', label=r'$t_{\rm{WD}}=70\, \rm{Myr}$')    

        self.ax2.plot(
          [np.nan], [np.nan], color=ct[2], ls=lst[2], lw=lw,
          marker='None', label=r'$t_{\rm{WD}}=100\, \rm{Myr}$')    

        self.ax2.legend(
          frameon=False, fontsize=fs, numpoints=1, ncol=1, labelspacing=.2,
          handlelength=1.5, handletextpad=.8, loc=3, bbox_to_anchor=(.2, 0.)) 


    def manage_output(self):
        if self.save_fig:
            fpath = './../OUTPUT_FILES/ANALYSES_FIGURES/Fig_sSNRL-tests.pdf'
            plt.savefig(fpath, format='pdf')
        if self.show_fig:
            plt.show() 
        plt.close()
            
    def make_plot(self):
        self.set_fig_frame()
        self.plot_models()
        self.make_legend()
        self.manage_output()

if __name__ == '__main__':
    Plot_Tests(show_fig=False, save_fig=True)
 
