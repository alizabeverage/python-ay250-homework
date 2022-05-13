import sys
import warnings
from copy import deepcopy
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy import constants, interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.table import Table, Column, hstack
from collections import OrderedDict

class Alf(object):
    def __init__(self, outfiles):
        self.outfiles = outfiles
        self.mcmc = np.load('{0}.npy'.format(self.outfiles))
        self.labels = np.array([
                  'chi2','velz','sigma','logage','zH',
                  'FeH', 'a', 'C', 'N', 'Na', 'Mg', 'Si',
                  'K', 'Ca', 'Ti','V', 'Cr', 'Mn', 'Co',
                  'Ni', 'Cu', 'Sr','Ba', 'Eu', 'Teff',
                  'IMF1', 'IMF2', 'logfy', 'sigma2', 'velz2',
                  'logm7g', 'hotteff', 'loghot', 'fy_logage',
                  'logtrans', 'logemline_h', 'logemline_oii',
                  'logemline_oiii', 'logemline_sii', 'logemline_ni',
                  'logemline_nii', 'jitter', 'IMF3', 'logsky', 'IMF4',
                  'h3', 'h4', 'ML_v','ML_i','ML_k','MW_v', 'MW_i','MW_k'
                  ])
        
        self.elements = np.array(['a', 'C', 'N', 'Na', 'Mg','Si',
                    'K', 'Ca', 'Ti','V', 'Cr','Mn', 
                    'Co', 'Ni', 'Cu', 'Sr','Ba','Eu'])

        self.results = self.get_abundances()
        for i,label in enumerate(self.labels):
            if label not in self.elements: 
                self.results[label] = np.around(np.percentile(self.mcmc[:,i],[16,50,84]),4)
           
        self.get_feh()
        self.results = Table(self.results)
        
#         self.spectra = self.get_spectrum()
        
        
        
        
    def get_spectrum(self):
        m = np.loadtxt('{0}.bestspec'.format(self.outfiles))
        data = {}
        data['wave'] = m[:,0]/(1.+self.results['velz'][1]*1e3/constants.c)
        data['m_flux'] = m[:,1] # Model spectrum, normalization applied
        data['d_flux'] = m[:,2] # Data spectrum
        data['snr'] = m[:,3]  # Including jitter and inflated errors
        data['unc'] = m[:,2]/m[:,3]
        data['poly'] = m[:,4] # Polynomial used to create m_flux
        data['residual'] = (m[:,1] - m[:,2])/m[:,1] * 1e2
        
        return data
    
    def get_feh(self):

        zh = np.where(self.labels == 'zH')
        feh = np.where(self.labels == 'FeH')
        total_met = self.mcmc[:,zh] + self.mcmc[:,feh]

        #Computing errors directly from the chains.
        self.results['feh'] = np.around(np.percentile(total_met,[16,50,84]),4)

    
    def get_abundances(self, get_distributions=False, s07=False, b14=False, m11=True):
        """
        Convert abundances from X/H to X/Fe.

        Correct the raw abundance values given
        by ALF.
        """
        # initialize
        abundances = OrderedDict()
        abundances['Type'] = ['16','50','84']
        distributions = OrderedDict()
        
        # Correction factros from Schiavon 2007, Table 6
        # NOTE: Forcing factors to be 0 for [Fe/H]=0.0,0.2
        lib_feh = [-1.6, -1.4, -1.2, -1.0, -0.8,
                   -0.6, -0.4, -0.2, 0.0, 0.2]
        lib_ofe = [0.6, 0.5, 0.5, 0.4, 0.3, 0.2,
                   0.2, 0.1, 0.0, 0.0]

        if s07:
            #Schiavon 2007
            lib_mgfe = [0.4, 0.4, 0.4, 0.4, 0.29,
                        0.20, 0.13, 0.08, 0.05, 0.04]
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.20,
                        0.12, 0.06, 0.02, 0.0, 0.0]
        elif b14:
            # Fitted from Bensby+ 2014
            lib_mgfe = [0.4 , 0.4, 0.4, 0.38, 0.37,
                        0.27, 0.21, 0.12, 0.05, 0.0]
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26,
                        0.17, 0.12, 0.06, 0.0, 0.0]
        elif m11 or (b14 is False and s07 is False):
            # Fitted to Milone+ 2011 HR MILES stars
            lib_mgfe = [0.4, 0.4, 0.4, 0.4, 0.34, 0.22,
                        0.14, 0.11, 0.05, 0.04]
            # from B14
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26,
                        0.17, 0.12, 0.06, 0.0, 0.0]

        # In ALF the oxygen abundance is used
        # a proxy for alpha abundance
        del_alfe = interpolate.interp1d(lib_feh, lib_ofe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
        del_mgfe = interpolate.interp1d(lib_feh, lib_mgfe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
        del_cafe = interpolate.interp1d(lib_feh, lib_cafe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')

        zh = np.where(self.labels == 'zH')
        al_corr = del_alfe(self.mcmc[:,zh])
        mg_corr = del_mgfe(self.mcmc[:,zh])
        ca_corr = del_cafe(self.mcmc[:,zh])

        # Assuming Ca~Ti~Si
        group1 = {'Ca', 'Ti', 'Si'}

        # These elements seem to show no net enhancemnt
        # at low metallicity
        group2 = {'C', 'Ca', 'N', 'Cr', 'Ni', 'Na'}

        # These elements we haven't yet quantified
        group3 = {'Ba', 'Eu', 'Sr', 'Cu', 'Co',
                  'K', 'V', 'Mn'}


        
    
        mgh = np.where(self.labels == 'Mg')
        feh = np.where(self.labels == 'FeH')
        zh = np.where(self.labels == 'zH')
        for i, element in enumerate(self.elements):
            xh = np.where(self.labels == element)
            xfe = (self.mcmc[:,xh] - self.mcmc[:,feh])
            if element=='Mg':
                xfe_vals = xfe + mg_corr
                xh_vals = self.mcmc[:,xh] + self.mcmc[:,zh] + mg_corr
                xmg_vals = self.mcmc[:,xh] - self.mcmc[:,mgh] + mg_corr - mg_corr
            elif element=='a':
                xfe_vals = xfe + al_corr
                xh_vals = self.mcmc[:,xh] + self.mcmc[:,zh] + al_corr
                xmg_vals = self.mcmc[:,xh] - self.mcmc[:,mgh] + al_corr - mg_corr
            elif element in group1:
                xfe_vals = xfe + ca_corr
                xh_vals = self.mcmc[:,xh] + self.mcmc[:,zh] + ca_corr
                xmg_vals = self.mcmc[:,xh] - self.mcmc[:,mgh] + ca_corr - mg_corr
            elif element in group2 or element in group3:
                xfe_vals = xfe
                xh_vals = self.mcmc[:,xh] + self.mcmc[:,zh]
                xmg_vals = self.mcmc[:,xh] - self.mcmc[:,mgh]

            # get distributions and save to ordered dict!
            xfe_dist = np.around(np.percentile(xfe_vals,[16,50,84]),4)
            xh_dist = np.around(np.percentile(xh_vals,[16,50,84]),4)
            xmg_dist = np.around(np.percentile(xmg_vals,[16,50,84]),4)
            
            abundances[element.lower()+'fe'] = xfe_dist
            abundances[element.lower()+'h'] = xh_dist
            abundances[element.lower()+'mg'] = xmg_dist
            
            if get_distributions:
                distributions[element.lower()+'fe'] = xfe_vals
                distributions[element.lower()+'h'] = xh_vals
                distributions[element.lower()+'mg'] = xmg_vals

                return distributions
            
        # remove 'mgmg' because that means nothing...
        abundances.pop('mgmg')
        return abundances

    def get_cls(self, distribution):
        distribution = np.sort(np.squeeze(distribution))
        lower,median,upper = np.percentile(distribution,[16,50,84])
 
        return {'cl16': np.round(lower,decimals=4), \
                'cl50':  np.round(median,decimals=4), \
                'cl84': np.round(upper,decimals=4)}

    def plot_corner(self, plot_labels, figure = None, color='C0'):
        sel_labels = [l in plot_labels for l in o.labels]
        corner_labels = o.labels[sel_labels]
        mcmc_sel = o.mcmc[:,sel_labels]

        figure = corner.corner(mcmc_sel, labels=corner_labels,
                       show_titles=True, title_kwargs={"fontsize": 12},
                       color=color,fig=figure)

        return figure


    
def make_results_tab(fnames,savename=None):
    '''
    fnames: include * in the filename, use glob to search for all files
    append all results to a table
    save if you want
    '''
    import glob
    files = glob.glob(fnames+'mcmc')
    files = [f[:-5] for f in files]

    for i,f in enumerate(files):
        o = Alf(f)
        if i == 0:
            # initialize dictionary
            od16 = OrderedDict({key+'16':[] for key in o.results.colnames[1:]})
            od50 = OrderedDict({key+'50':[] for key in o.results.colnames[1:]})
            od84 = OrderedDict({key+'84':[] for key in o.results.colnames[1:]})
            od = {'file':[],**od16, **od50, **od84}

        od['file'].append(f.split('/')[-1])
        [od[label+'16'].append(o.results[label][0]) for label in o.results.colnames[1:]]
        [od[label+'50'].append(o.results[label][1]) for label in o.results.colnames[1:]]
        [od[label+'84'].append(o.results[label][2]) for label in o.results.colnames[1:]]

        del o
    
    results = Table(od)
    if savename is not None:
        results.write(savename+'.fits',overwrite=True)
        
    return results

def correct_abundance(zh,xh,element):
    if element not in ['cah', 'tih', 'sih','mgh','ah']:
        return xh+zh
    
    lib_feh = [-1.6, -1.4, -1.2, -1.0, -0.8,
                   -0.6, -0.4, -0.2, 0.0, 0.2]
    lib_ofe = [0.6, 0.5, 0.5, 0.4, 0.3, 0.2,
                   0.2, 0.1, 0.0, 0.0]
    lib_mgfe = [0.4, 0.4, 0.4, 0.4, 0.34, 0.22,
                   0.14, 0.11, 0.05, 0.04]
    lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26,
                   0.17, 0.12, 0.06, 0.0, 0.0]
    
    del_alfe = interpolate.interp1d(lib_feh, lib_ofe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
    del_mgfe = interpolate.interp1d(lib_feh, lib_mgfe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
    del_cafe = interpolate.interp1d(lib_feh, lib_cafe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
    
    
    al_corr = del_alfe(zh)
    mg_corr = del_mgfe(zh)
    ca_corr = del_cafe(zh)
    
    if element=='mgh':
        return xh + zh + mg_corr

    elif element=='ah':
        return xh + zh + al_corr
    
    elif element in ['cah', 'tih', 'sih']:
        return xh + zh + ca_corr
    
    print("something's not write")

    
    
    
        
    
    
    
