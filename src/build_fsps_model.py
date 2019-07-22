#!/usr/bin/env python
import sys, os, time
import numpy as np
import pandas as pd
from astropy import units as u

class Build_Fsps(object):
    """
    Description:
    ------------
    TBW.
    Parameters:
    -----------
    _inputs : ~instance
        Instance of the Input_Parameters class defined in input_params.py.
    _tau : ~float
        Timescale.
     
    Outputs:
    --------
    TBW
    """     
    def __init__(self, _inputs):
        self._inputs = _inputs
        
        self.D = {}
        self.D['Dcd_fine'] = np.arange(-1.1, 1.00001, 0.01)
        self.D['t_ons'] = self._inputs.t_onset.to(u.Gyr).value
        self.D['t_bre'] = self._inputs.t_cutoff.to(u.Gyr).value

        #Improve this.
        #abs_dir = os.path.split(os.path.abspath(os.path.realpath(sys.argv[0])))[0]
        try:
            synpop_dir = os.path.join(os.environ['PATH_ssnarl'] + 'fsps_files/' + self._inputs.fsps_path)
        except:
            synpop_dir = os.path.join('./../fsps_files/' + self._inputs.fsps_path)

        #Get SSP data and compute the theoretical color with respect to the RS.
        df = pd.read_csv(synpop_dir + 'SSP.dat', header=0, escapechar='#')
        logage_SSP = df[' log_age'].values
        mag_2_SSP = df[self._inputs.f2].values
        mag_1_SSP = df[self._inputs.f1].values
        
        #Retrieve RS color.
        RS_condition = (logage_SSP == 10.0)
        RS_color = mag_2_SSP[RS_condition] - mag_1_SSP[RS_condition]

        self.D['logage'] = logage_SSP
        for i, tau in enumerate(self._inputs.tau_list):
            TS = str(tau.to(u.yr).value / 1.e9)
         
            model = pd.read_csv(synpop_dir + 'tau-' + TS + '.dat', header=0)
            self.D['tau_' + TS] = tau.to(u.Gyr).value
            self.D['mag2_' + TS] = model[self._inputs.f2].values
            self.D['mag1_' + TS] = model[self._inputs.f1].values
            self.D['mag0_' + TS] = model[self._inputs.f0].values
            self.D['age_' + TS] = 10.**(model['# log_age'].values) / 1.e9 #Converted to Gyr.
            self.D['int_mass_' + TS] = model['integrated_formed_mass'].values
            self.D['int_stellar_mass_' + TS] = model['integrated_stellar_mass'].values
            self.D['Dcolor_' + TS] = (
              self.D['mag2_' + TS] - self.D['mag1_' + TS] - RS_color) 

            #Get analytical normalization for the SFH.
            if self._inputs.sfh_type == 'exponential':
                self.D['sfr_norm_' + TS] = (
                  -1. / (self.D['tau_' + TS] * (np.exp(-self.D['age_' + TS][-1] / self.D['tau_' + TS])
                  - np.exp(-self.D['age_' + TS][0] / self.D['tau_' + TS]))))     
            elif self._inputs.sfh_type == 'delayed-exponential':
                self.D['sfr_norm_' + TS] = (
                  1. / (((-self.D['tau_' + TS] * self.D['age_' + TS][-1] - self.D['tau_' + TS]**2.)
                  * np.exp(- self.D['age_' + TS][-1] / self.D['tau_' + TS])) -
                  ((-self.D['tau_' + TS] * self.D['age_' + TS][0] - self.D['tau_' + TS]**2.) * np.exp(- self.D['age_' + TS][0] / self.D['tau_' + TS]))))
