import numpy as np
import h5py
from numpy.polynomial import chebyshev
import pandas as pd
from scipy import interpolate

#calculation of the electron capture/charge exchnage cross section between proton and hydrogen molecules, assuming beam-target approach (T_H2 = 0).

#measured data: BARNETT, Clarence F., et al. Atomic data for fusion. Volume 1: Collisions of H, H2, He and Li atoms and ions with atoms and molecules. NASA STI/Recon Technical Report N, 1990, 91: 13238.
#measued data 4s: HUGHES, R. H.; DAWSON, H. R.; DOUGHTY, B. M. Electron Capture into the 4 s State of H by Fast H+ Impact on Gases. Physical Review, 1967, 164.1: 166.
#theoretical data: PLOWMAN, Corey T., et al. Effective one-electron approach to proton collisions with molecular hydrogen. The European Physical Journal D, 2022, 76.2: 31.(unpublished states from the personal communication)

#Measured data are based on the chebyshev fit from Bernett 1990
#n=1: 2.6e-3 - 4e3 keV, Barnett
#n=2: 8.2e-2 - 46.984 keV, Barnett, Lyman excitation, 46.984 - 71.660 keV, Barnett, 2s+2p, 71.660 - 2000 keV Plowman
#n=3: 2.2 - 110.909 keV Barnett, 3s+3p+3d, 110.909 - 2000 keV Plowman
#n=4: 4s+4p+4d (4f not considered), 4s: 5-115.860 keV, Hughes 1967, 115.860-2000 keV Plowman, 4p, 4d: 5-2000 keV Plowman
#n=5,6: 5-2000 keV Plowman (5,6f not considered)


class calculate_H_H2_FIDASIM:
    #nenergy: number of energy points, Emin: min energy of calculated range in keV, Emax: max of calculated range in keV, n_end: principal quantum numbers of the end state of new atom, input_files_path: file with cross-sections parameters, file_to_save: h5 file where the parameters will be saved, cross section in cm2
    
    def __init__(self, nenergy, Emin, Emax, n_end = np.arange(1,7), input_files_path = '', save_file = True, file_to_save = 'atomic_tables.h5'):
        emin_cx, emax_cx = 5, 2000 #energy range for h+h2 cx cross-section available for all n=1-6
        if Emin<emin_cx:
            print('too low Emin for CX H+H2, minimum is {} keV'.format(emin_cx))
        else:
            pass
        if Emax>emax_cx:
            print('too high Emax for CX H+H2, maximum is {} keV'.format(emax_cx))
        else:
            pass
        
        #name of folders with cx H-H2 cross section based on measuremets
        Emin_log, Emax_log = np.log10(Emin), np.log10(Emax)
        dE_log = (Emax_log-Emin_log)/(nenergy-1)
        E_log = np.arange(Emin_log, Emax_log+dE_log, dE_log)
        sigma = np.zeros((len(n_end), len(E_log)))
        
        for n in n_end:
            sigma[n-1, :] = self.cx_sigma_calc(E_log, n, input_files_path)
        
        sigma = np.nan_to_num(sigma)
            
        if save_file:
            self.save_file(Emin, Emax, dE_log, E_log, sigma, file_to_save)
            print('cx H+H2 saved into {}'.format(file_to_save))
        else:
            self.sigma = sigma
            self.E_log = E_log
        
    
    def Chebyshev_Barnett(self, f_s, E): #E in eV,  Emin and Emax are ranges of available cross sections, E is required energy range
        Emin_s, Emax_s = f_s['E_min'][()], f_s['E_max'][()] #limits in eV
        cheb_coef = np.array([f_s[name][()] for name in f_s.keys() if "A" in name])
            
        E_new = np.array([ee for ee in E if Emin_s<=ee<=Emax_s])
        if E_new[0] > E[0]:
            print('Emin is lower than defined range for p-H2 CX collisions')
        else: pass
            
        if E_new[-1] > E[-1]:
            print('Emax is higher than defined range for p-H2 CX collisions')
        else: pass
            
        x = ((np.log(E_new)-np.log(Emin_s))-(np.log(Emax_s)-np.log(E_new)))/(np.log(Emax_s)-np.log(Emin_s))
        y = chebyshev.chebval(x, cheb_coef)-cheb_coef[0]/2
        
        cross_sec = np.zeros_like(E)
        cross_sec[np.isin(E,E_new)] = np.exp(y)
        
        
        #return sigma in cm2, E in keV
        return E_new/1000.,cross_sec
    
    def spline(self, ee, ss, E):
        
        p = interpolate.interp1d(np.log10(ee), np.log10(ss), bounds_error=False)
        return 10**p(np.log10(E))
    
    def cx_sigma_calc(self, E_log, n_end, input_files_path):
        print('calculate cx H+H2 for n_end {}'.format(n_end))
        Emin, Emax = 10**(E_log[0]), 10**(E_log[-1])
    
        measured = {1:['h_h2_n1s.h5'], 
                    2:['h_h2_n2s.h5', 'h_h2_n2p.h5', 'h_h2_n2.h5'], 
                    3:['h_h2_n3s.h5', 'h_h2_n3p.h5', 'h_h2_n3d.h5'], 
                    4:['4s_measured_Hughes.csv'], 
                    5:[np.nan], 
                    6:[np.nan]}
        
        #name of folders with cx p-H2 cross section based on theory
        Plowman = {1:'ec1.dat',  
                   2:'ec2.dat', 
                   3:'ec3.dat', 
                   4:'ec4.dat',
                   5:'ec5.dat', 
                   6:'ec6.dat' }     
        
        #calculating p-H2 cross-section for end state n=1
        if n_end == 1:
            print('defined energy range of H+H2 for n_end=1: 2.6e-3 - 4e3 keV')
            f_name = measured[n_end][0]
            f_s = h5py.File(input_files_path+'measured_H_H2/'+f_name,'r+')
            E, sigma = self.Chebyshev_Barnett(f_s, 1000*10**E_log)
            f_s.close()
            
        #calculating p-H2 cross-section for end state n=2  
        elif n_end == 2:
            print('defined energy range of H+H2 for n_end=2: 8.2e-2 - 2e3 keV')
            #based on p_h2_n2 up to E_b1, sum of p_h2_n2s and p_h2_n2p between E_b1, E_b2, based on Plowmann above E_b2
            E_b1, E_b2 = 46.984, 71.660 #keV
            f_names = measured[n_end]
            f_names_th = Plowman[n_end]

            if Emin < E_b1:
                print('calculate sigma low')
                #reading parameters of Chebyshev polynoms of measured data fits
                f_2 = h5py.File(input_files_path+'measured_H_H2/'+f_names[2],'r+')

                E_low = 10**np.array([E for E in E_log if 10**E < E_b1])
                sigma_low = self.Chebyshev_Barnett(f_2, 1000*E_low)
                f_2.close()
            else:
                sigma_low = np.array([]), np.array([])

            if Emin <= E_b2 and Emax >= E_b1:
                print('calculate sigma mid')
                #reading parameters of Chebyshev polynoms of measured data fits
                f_2s = h5py.File(input_files_path+'measured_H_H2/'+f_names[0],'r+') 
                f_2p = h5py.File(input_files_path+'measured_H_H2/'+f_names[1],'r+')

                E_mid = 10**np.array([E for E in E_log if (10**E >= E_b1 and 10**E <= E_b2)])
                sigma_2s = self.Chebyshev_Barnett(f_2s, 1000*E_mid)
                sigma_2p = self.Chebyshev_Barnett(f_2p, 1000*E_mid)
                sigma_mid = E_mid, sigma_2s[1] + sigma_2p[1]
                f_2s.close(), f_2p.close()
            else:
                sigma_mid = np.array([]), np.array([])

            if Emax > E_b2:
                print('calculate sigma high')
                #reading theoretical data
                n2 = pd.read_csv(input_files_path+'calculated_H_H2/'+f_names_th, names = ['E', 's', 'p'], sep = '    ')
                E_high = 10**np.array([E for E in E_log if 10**E > E_b2])
                sigma_high = E_high, self.spline(n2.E, (n2.s+n2.p)*10**(-16), E_high)
            else:
                sigma_high = np.array([]), np.array([])

            #E = np.concatenate((sigma_low[0], sigma_mid[0], sigma_high[0]))
            sigma = np.concatenate((sigma_low[1], sigma_mid[1], sigma_high[1]))

        #calculating p-H2 cross-section for end state n=3    
        elif n_end == 3:
            print('defined energy range of H+H2 for n_end=3: 2.2 - 2e3 keV')
            #based on measurement up to E_b, based on Plowmann above E_b
            E_b = 110.909 #keV
            f_names = measured[n_end]
            f_names_th = Plowman[n_end]

            if Emin <= E_b:
                print('calculate sigma low')
                #reading parameters of Chebyshev polynoms of measured data fits
                f_3s = h5py.File(input_files_path+'measured_H_H2/'+f_names[0],'r+') 
                f_3p = h5py.File(input_files_path+'measured_H_H2/'+f_names[1],'r+')
                f_3d = h5py.File(input_files_path+'measured_H_H2/'+f_names[2],'r+')

                E_low = 10**np.array([E for E in E_log if (10**E <= E_b)])
                sigma_3s = self.Chebyshev_Barnett(f_3s, 1000*E_low)
                sigma_3p = self.Chebyshev_Barnett(f_3p, 1000*E_low)
                sigma_3d = self.Chebyshev_Barnett(f_3d, 1000*E_low)

                sigma_low = E_low, sigma_3s[1] + sigma_3p[1] + sigma_3d[1]
                f_3s.close(), f_3p.close(), f_3d.close()
            else:
                sigma_low = np.array([]), np.array([])

            if Emax > E_b:
                print('calculate sigma high')
                #reading theoretical data
                n3 = pd.read_csv(input_files_path+'calculated_H_H2/'+f_names_th, names = ['E', 's', 'p', 'd'], sep = '    ')
                E_high = 10**np.array([E for E in E_log if 10**E > E_b])
                sigma_high = E_high, self.spline(n3.E, (n3.s+n3.p+n3.d)*10**(-16), E_high)
            else:
                sigma_high = np.array([]), np.array([])

            #E = np.concatenate((sigma_low[0], sigma_high[0]))
            sigma = np.concatenate((sigma_low[1], sigma_high[1]))

        #calculating p-H2 cross-section for end state n=4
        elif n_end == 4:
            print('defined energy range of H+H2 for n_end=4: 5 - 2e3 keV')
            E_b = 115.860 #keV
            f_names = measured[n_end]
            f_names_th = Plowman[n_end]

            n4 = pd.read_csv(input_files_path+'calculated_H_H2/'+f_names_th, names = ['E', 's', 'p', 'd'], sep = '    ')
            sigma_4pd = self.spline(n4.E, (n4.p+n4.d)*10**(-16), 10**E_log)

            if Emin <= E_b:
                n4s_m = pd.read_csv(input_files_path+'measured_H_H2/'+f_names[0], names = ['E', 's'], sep = ',')
                E_low = 10**np.array([E for E in E_log if 10**E <= E_b])
                sigma_4s_low = E_low, self.spline(n4s_m.E, n4s_m.s*10**(-19), E_low)
            else:
                sigma_4s_low = np.array([]), np.array([])

            if Emax > E_b:
                print('calculate sigma high')
                #reading theoretical data
                n4s = pd.read_csv(input_files_path+'calculated_H_H2/'+f_names_th, names = ['E', 's', 'p', 'd'], sep = '    ')
                E_high = 10**np.array([E for E in E_log if 10**E > E_b])
                sigma_4s_high = E_high, self.spline(n4.E, (n4.s)*10**(-16), E_high)
            else:
                sigma_4s_high = np.array([]), np.array([])

            #E = 10**E_log
            sigma = np.concatenate((sigma_4s_low[1], sigma_4s_high[1])) + sigma_4pd

        #calculating p-H2 cross-section for end state n=5   
        elif n_end == 5:
            print('defined energy range of H+H2 for n_end=5: 5 - 2e3 keV')
            f_names_th = Plowman[n_end]
            n5 = pd.read_csv(input_files_path+'calculated_H_H2/'+f_names_th, names = ['E', 's', 'p', 'd'], sep = '    ')
            #E = 10**E_log
            sigma = self.spline(n5.E, (n5.s+n5.p+n5.d)*10**(-16), 10**E_log)

        #calculating p-H2 cross-section for end state n=6    
        elif n_end == 6:
            print('defined energy range of H+H2 for n_end=6: 5 - 2e3 keV')
            f_names_th = Plowman[n_end]
            n6 = pd.read_csv(input_files_path+'calculated_H_H2/'+f_names_th, names = ['E', 's', 'p', 'd'], sep = '    ')
            #E = 10**E_log
            sigma = self.spline(n6.E, (n6.s+n6.p+n6.d)*10**(-16), 10**E_log)
        else:
            print('incorrect end state')
        return sigma
    
    
    def save_file(self, Emin, Emax, dE_log, E_log, sigma, file_to_save):
        with h5py.File(file_to_save, 'a') as f:

            p_H2 = f.create_group('/cross/H_H2')
            p_H2.attrs['description'] = 'Cross sections for Hydrogen-Hydrogen molecule interactions'

            crs_d = f.create_dataset('/cross/H_H2/cx', data = sigma.T)
            crs_d.attrs['description'] = 'm resolved charge exchange cross sections: cx(m,energy)'
            crs_d.attrs['reaction'] = 'H(+) + H2 -> H(m) + ...'
            crs_d.attrs['units'] = 'cm^2'

            dloge_d = f.create_dataset('/cross/H_H2/dlogE', data = np.round(dE_log, 7))
            dloge_d.attrs['description'] = 'Energy spacing in log-10'
            dloge_d.attrs['units'] = 'log10(keV/amu)'

            emax_d = f.create_dataset('/cross/H_H2/emax', data = Emax)
            emax_d.attrs['description'] = 'Maximum energy'
            emax_d.attrs['units'] = 'keV/amu'

            emin_d = f.create_dataset('/cross/H_H2/emin', data = Emin)
            emin_d.attrs['description'] = 'Minimim energy'
            emin_d.attrs['units'] = 'keV/amu'

            nenergy_d = f.create_dataset('/cross/H_H2/nenergy', data = len(E_log))
            nenergy_d.attrs['description'] = 'Number of nucleon energy values'

            energy_d = f.create_dataset('/cross/H_H2/energy', data = 10**E_log)
            energy_d.attrs['description'] = 'Nucleon energy values'
            energy_d.attrs['units'] = 'keV/amu'

            m_max_d = f.create_dataset('/cross/H_H2/m_max', data = sigma.shape[0])
            m_max_d.attrs['description'] ='Number of final energy levels'
        
