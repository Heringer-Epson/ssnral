import numpy as np
from astropy import units as u

class Generic_Pars(object):
    """
    Description:
    ------------
    Define a set of input parameters to use to make the Dcolor vs SN rate plot
    in the class below. Detailed description of the parameters to follow.

    Parameters:
    -----------
    TBW.
    """  
    def __init__(
      self, sfh_type='exponential', imf_type='Kroupa', Z='0.0190',
      t_onset=1.e8*u.yr, t_cutoff=1.e9*u.yr, fbhb=0.0, dust='0',
      spec_lib='BASEL', isoc_lib='PADOVA', f2='g', f1='r', f0='r'):

        self.sfh_type = sfh_type
        self.imf_type = imf_type
        self.Z = Z
        self.t_onset = t_onset
        
        self.t_cutoff = t_cutoff
        self.fbhb = fbhb
        self.dust = dust
        self.spec_lib = spec_lib
        self.isoc_lib = isoc_lib
        self.f2, self.f1, self.f0 = f2, f1, f0

        #Do not change the following.
        self.tau_list = np.array(
          [1., 1.5, 2., 3., 4., 5., 7., 10.]) * 1.e9 * u.yr
        self.subdir_fullpath = './'
                
        self.fsps_path = (
          self.imf_type + '_' + self.sfh_type + '_' + self.Z + '_'
          + str(self.fbhb) + '_' + self.dust + '_' + self.spec_lib + '_'
          + self.isoc_lib + '/')
