import os, shutil
import numpy as np
from astropy import units as u

class Generic_Pars(object):
    """
    Description:
    ------------
    Define a set of input parameters to use to make the Dcolor vs SN rate plot
    in the class below. This is intended to replicate the ./../src/ code
    input_params.py, but only containing the relevant quantities for this plot.

    Parameters:
    -----------
    As described in ./../src/input_params.py
    """  
    def __init__(
      self, sfh_type, imf_type, Z, t_onset, 
      t_cutoff=1.e9*u.yr, fbhb=0.0, spec_lib='BASEL', isoc_lib='PADOVA'):

        self.sfh_type = sfh_type
        self.imf_type = imf_type
        self.Z = Z
        self.t_onset = t_onset
        
        self.t_cutoff = t_cutoff
        self.spec_lib = spec_lib
        self.isoc_lib = isoc_lib
        self.fbhb = fbhb

        #Do not change the following.
        self.filter_1 = 'r'
        self.filter_2 = 'g'
        self.tau_list = np.array(
          [1., 1.5, 2., 3., 4., 5., 7., 10.]) * 1.e9 * u.yr
        self.subdir_fullpath = './'
                
        self.fsps_path = (
          self.imf_type + '_' + self.sfh_type + '_' + self.Z + '_'
          + str(self.fbhb) + '_' + self.spec_lib + '_' + self.isoc_lib + '/')
