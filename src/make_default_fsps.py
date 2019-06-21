#!/usr/bin/env python

import os
import time
import fsps
import numpy as np
from astropy import units as u

sfh_type2sfh_key = {'exponential': 1, 'delayed-exponential': 4}
imf_type2imf_key = {'Salpeter': 0, 'Chabrier': 1, 'Kroupa': 2}
Z2Z_key_PADOVA = {
  '0.0002': 1, '0.0003': 2, '0.0004': 3, '0.0005': 4, '0.0006': 5, '0.0008': 6,
  '0.0010': 7, '0.0012': 8, '0.0016': 9, '0.0020': 10, '0.0025': 11, '0.0031': 12,
  '0.0039': 13, '0.0049': 14, '0.0061': 15, '0.0077': 16, '0.0096': 17,
  '0.0120': 18, '0.0150': 19, '0.0190': 20, '0.0240': 21, '0.0300': 22}
Z2Z_key_BASTI = {
  '0.0003': 1, '0.0006': 2, '0.0010': 3, '0.0020': 4, '0.0040': 5, '0.0080': 6,
  '0.0100': 7, '0.0200': 8, '0.0300': 9, '0.0400': 10}
Z2Z_key_GENEVA = {
  '0.0010': 1, '0.0040': 2, '0.0080': 3, '0.0200': 4, '0.0400': 5}
Z2Z_key_PARSEC = {
  '0.0001': 1, '0.0002': 2, '0.0005': 3, '0.0010': 4, '0.0020': 5, '0.0040': 6,
  '0.0060': 7, '0.0080': 8, '0.0100': 9, '0.0140': 10, '0.0170': 11, '0.0200': 12,
  '0.0300': 13, '0.0400': 14, '0.0600': 15}
Z2Z_key_MIST = {'0.0190': 10}

class Make_FSPS(object):
    """
    Description:
    ------------
    This is a stand alone routine to generate default FSPS files which can be
    used in case the user does have fsps and python-fsps installed.

    Parameters:
    -----------
    imf : ~str
        Initial mass function. 'Salpeter', 'Chabrier' or 'Kroupa'.
        Default in H17 is 'Chabrier'.
    sfh : ~str
        Star formation history. 'exponential' or 'delayed-exponential'
        Default in H17 is 'exponential'.
    Z : ~float
        Metallicity of the stellar population.
        Default in H17 is 0.0190 (solar).
    dust : ~int
        Dust attenuation around old stars, default is 0 (I checked this).
    spec_lib : ~str
        'BASEL' or 'MILES'.
        Default in H17 is 'BASEL'.
        IMPORTANT: This flag needs to be set in sps_vars.f90 in the src/
        directory of where FSPS is installed and FSPS needs to be re-compiled.
        It's function here is to then allocate the generated files in the
        correct sub-directories.
    isoc_lib : ~str
        'PADOVA', 'PARSEC', 'BASTI', GENEVA'.
        Default in H17 is 'PADOVA'.
        IMPORTANT: This flag needs to be set in sps_vars.f90 in the src/
        directory of where FSPS is installed and FSPS needs to be re-compiled.
        It's function here is to then allocate the generated files in the
        correct sub-directories. Notice this matters for setting the proper
        metallicity key.      
     
    Outputs:
    --------
    ./../fsps_files/fname.dat
        where X is the chosen SFH ('exponential' or 'delayed-exponential')
        and Y in the timescale.

    Notes:
    ------
    For new calculations that use different spectral or isochrone libraries,
    FSPS needs to be recompiled and Python-FSPS needs to be reinstalled.

    References:
    -----------
    Heringer+ 2017 (H17): http://adsabs.harvard.edu/abs/2017ApJ...834...15H
    """  

    def __init__(self, imf, sfh, Z, fbhb, dust, spec_lib, isoc_lib):
        print '\n\n>RUNNING FSPS... (requires fsps and python-fsps installed!)\n'
        self.filter_1 = 'r'
        self.filter_2 = 'g'
        self.imf_type = imf
        self.sfh_type = sfh
        self.Z = Z
        self.fbhb = fbhb
        self.dust = dust
        self.spec_lib = spec_lib
        self.isoc_lib = isoc_lib
                
        _tau_list = [1., 1.5, 2., 3., 4., 5., 7., 10.]
        self.tau_list = [tau * 1.e9 * u.yr for tau in _tau_list]    
        self.directory = None
          
        self.make_directory()
        self.run_fsps() 

    def make_directory(self):
        self.directory = (
          './../fsps_files/' + self.imf_type + '_' + self.sfh_type
          + '_' + self.Z  + '_' + str(self.fbhb) + '_' + str(self.dust) + '_'
          + self.spec_lib + '_' + self.isoc_lib + '/')
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)             
        
    def make_output(self, sp, fname):
        """Given a FSPS object (sp), write an output file."""
        header = ('log_age,instantaneous_sfr,integrated_stellar_mass,'\
                  + 'integrated_formed_mass,u,g,r,i,z')

        #Data needs to be transposed, otherwise arrays are written as
        #lines rather than columns.
        mag_u, mag_g, mag_r, mag_i, mag_z = zip(*sp.get_mags(
          bands=['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']))
        out_data = (sp.log_age, sp.sfr, sp.stellar_mass, sp.formed_mass,
          mag_u, mag_g, mag_r, mag_i, mag_z)
        np.savetxt(self.directory + fname, np.transpose(out_data), header=header,
          delimiter=',')          

    def run_fsps(self):
        
        imf_key = imf_type2imf_key[self.imf_type]
        sfh_key = sfh_type2sfh_key[self.sfh_type]
        if self.isoc_lib == 'PADOVA':
            Z_key = Z2Z_key_PADOVA[self.Z]
        elif self.isoc_lib == 'BASTI':
            Z_key = Z2Z_key_BASTI[self.Z]
        elif self.isoc_lib == 'GENEVA':
            Z_key = Z2Z_key_GENEVA[self.Z]
        elif self.isoc_lib == 'PARSEC':
            Z_key = Z2Z_key_PARSEC[self.Z]
        elif self.isoc_lib == 'MIST':
            Z_key = Z2Z_key_MIST[self.Z]

        #Make Simple Stellar Population (SSP).
        print '  *Calculating a SSP.'
        sp = fsps.StellarPopulation(
          compute_vega_mags=False, zcontinuous=0, zmet=Z_key,
          add_agb_dust_model=False, add_dust_emission=False,
          add_stellar_remnants=True, fbhb=self.fbhb, pagb=0., zred=0.,
          imf_type=imf_key, sfh=0, dust_type=self.dust, tage=14.96)        

        #Write output file.
        fname = 'SSP.dat'
        self.make_output(sp, fname)
    
        for tau in self.tau_list:
            tau = tau.to(u.yr).value / 1.e9
            print '    Tau=' + str(tau) + ' Gyr.'
            #Call FSPS to compute synthetic stellar populations 
            sp = fsps.StellarPopulation(
              compute_vega_mags=False, zcontinuous=0, zmet=Z_key,
              add_agb_dust_model=False, add_dust_emission=False,
              add_stellar_remnants=True, fbhb=self.fbhb, pagb=0., zred=0., imf_type=imf_key,
              sfh=sfh_key, tau=tau, const=0., fburst=0., dust_type=self.dust, tage=14.96)

            #Write output files.
            fname = 'tau-' + str(tau) + '.dat'
            self.make_output(sp, fname)            

        #Make record file.
        with open(self.directory + 'fsps_info.txt', 'w') as out:
            out.write('Files created on: ' + time.strftime('%d/%m/%Y') + '\n')
            out.write('FSPS version: ' + str(fsps.__version__) + '\n')
            out.write('Isochrones: ' + sp.libraries[0] + '\n')
            out.write('Spectral library: ' + sp.libraries[1] + '\n')
            out.write('IMF: ' + self.imf_type + '\n')
            out.write('SFH: ' + self.sfh_type + '\n')
            out.write('Metallicity: ' + str(self.Z))
        
if __name__ == '__main__':
    
    #Default.
    Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0190', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    
    #Dust.
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0190', fbhb=0.0, dust=1, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0190', fbhb=0.0, dust=2, spec_lib='BASEL', isoc_lib='PADOVA')
    
    #Delayed-exponential SFH.
    #Make_FSPS(imf='Kroupa', sfh='delayed-exponential', Z='0.0190', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    
    #Metallicity.
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0096', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0150', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0300', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    
    #Fraction of blue horizontal branch stars.
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0190', fbhb=0.1, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0190', fbhb=0.2, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
        
    #IMF.
    #Make_FSPS(imf='Chabrier', sfh='exponential', Z='0.0190', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Salpeter', sfh='exponential', Z='0.0190', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    
    #Additional tracks with Chabrier and Salpeter IMFs for other than solar Z.
    #Make_FSPS(imf='Chabrier', sfh='exponential', Z='0.0096', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Chabrier', sfh='exponential', Z='0.0300', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Salpeter', sfh='exponential', Z='0.0096', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')
    #Make_FSPS(imf='Salpeter', sfh='exponential', Z='0.0300', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='PADOVA')

    #=-=-=-=-=-=-=-=-=-=- RECOMPILE FSPS BEFORE RUNNING THESE -=-=-=-=-=-=-=-=-
    #Instructions:
    #In $SPS_HOME/src/:
    #Change the sps_vars.f90 accordingly (define BASEL = 0/1, define PADOVA=0/1)
    #>>make clean >>make
    #In $PY_FSPS_DIR:
    #>>python setup.py install
    #Then run this code, noting that the inputs for the spectral and isochrone
    #libraries must be consistent.
    
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0190', fbhb=0.0, dust=0, spec_lib='MILES', isoc_lib='PADOVA')    
    #Make_FSPS(imf='Kroupa', sfh='exponential', Z='0.0190', fbhb=0.0, dust=0, spec_lib='BASEL', isoc_lib='MIST')
