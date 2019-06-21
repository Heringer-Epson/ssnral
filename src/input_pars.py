import numpy as np
from astropy import units as u

class Pars(object):
    """
    Description:
    ------------
    Define a set of input parameters to use to make the Dcolor vs SN rate plot
    in the class below. Detailed description of the parameters to follow.

    Parameters:
    -----------
    imf_type : ~str
        Choice of initial mass function for the fsps simulations.
        Accepts: 'Salpeter', 'Chabrier' or 'Kroupa'. Default is 'Kroupa'.
    sfh_type : ~str
        Choice of star formation history for the fsps simulations.
        Accepts: 'exponential' (sfh propto exp(-t/tau)) or
        'delayed-exponential' (sfh propto t*exp(-t/tau).
        Default is 'exponential'.
    Z : ~str
        Choice of metallicity for the fsps simulations. Default is '0.0190'.
        Different choices of metallicity will require FSPS to be available.
    t_onset : ~astropy float (unit of time)
        Sets the time past which SN may occur. Determined by the lifetime of
        stars that may contribute to the SN rate. Default is .1*u.Gyr.
    t_cutoff : ~astropy float (unit of time)
        Sets the time at which the slope of the DTD may change. Determined by
        theoretical models. Default is 1.*u.Gyr.
    fbhb : float
        Fraction of blue horizontal branch stars. Options are 0.0 and 0.2, 
        depending on the availability of FSPS files. Default is 0.0.
    dust : int
        Dust treatment in FSPS. Defaul is 0.
    spec_lib : ~str
        Which spectral library to be used for FSPS simulations. Options are
        'BASEL' or 'MILES', depending on the availability of FSPS files.
        Default is 'BASEL'.
    isoc_lib_lib : ~str
        Which isochrone library to be used for FSPS simulations. Options are
        'PADOVA' or 'MIST', depending on the availability of FSPS files.
        Default is 'PADOVA'.      
    f2 : ~str
        filter_1 and filter_2 determine the color to be used as
        (filter_2 - filter_1). A list of available filters can be shown by
        calling fsps.list_filters() under run_fsps.py. Only g and r are
        currently supported. Default is 'g'.
    f1 : ~str
        Same as above. Default is 'r'.
    f0 : ~str
        The band used to for the luminosity unit. Default is 'r'.
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
