#!/usr/bin/env python

import numpy as np
from astropy import units as u

from input_pars import Pars
from build_fsps_model import Build_Fsps
from Dcolor2sSNRL_gen import Generate_Curve
from SN_rate import Model_Rates

class Produce_Function(object):
    """
    Description:
    ------------
    Provides an example of how to use the ssnral package.

    Parameters:
    -----------
    A : ~float
        DTD scale factor.
    s1 : ~float
        DTD slope at early times.
    s2 : ~float
        DTD slope at late times.        
    t_onset : ~astropy float (unit of time)
        Sets the time past which SN may occur. Determined by the lifetime of
        stars that may contribute to the SN rate.
    t_cutoff : ~astropy float (unit of time)
        Sets the time at which the slope of the DTD may change. Determined by
        theoretical models.       
    Z : ~float
        Metallicity. Used when computing FSPS files. Only a few options have
        synthetic stellar populations pre-computed here.
    imf_type : ~str
        Choice of initial mass function for the fsps simulations.
        Accepts: 'Salpeter', 'Chabrier' or 'Kroupa'.
    sfh_type : ~str
        Choice of star formation history for the fsps simulations.
        Accepts: 'exponential' (sfh propto exp(-t/tau)) or
                 'delayed-exponential' (sfh propto t*exp(-t/tau).

    Outputs:
    --------
    ./../OUTPUT_FILES/example_call.csv
    References:
    -----------
    Heringer+ 2017: http://adsabs.harvard.edu/abs/2017ApJ...834...15H
    """         
    def __init__(self, A, s1, s2, t_onset, t_cutoff, Z, sfh, imf):

        #Create a class containing input parameters.
        _inputs = Pars(
          sfh_type=sfh, imf_type=imf, Z=Z, t_onset=t_onset, t_cutoff=t_cutoff)
        
        #Computes quantities which require FSPS photometry, such as Delta(g-r).
        _D = Build_Fsps(_inputs).D 
        
        #Generate curves of the specific SN Ia rate as a function of color.
        Sgen = Generate_Curve(_inputs, _D, s1, s2)
        
        x,y = Sgen.Dcd_fine, np.log10(Sgen.sSNRL_fine * A)
        
        #Save output.
        fpath = './../OUTPUT_FILES/example_call.csv'
        np.savetxt(fpath, np.transpose((x,y)), fmt='%.4e', delimiter=',', header='D(g-r),sSNRL')
        print 'File saved at: ', fpath

if __name__ == '__main__':
    Produce_Function(
      A=1.e-12, s1=-1., s2=-2., t_onset=.1*u.Gyr, t_cutoff=.1*u.Gyr,
      Z='0.0190', sfh='exponential', imf='Kroupa')
