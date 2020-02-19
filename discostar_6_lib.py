import numpy as np
from scipy.interpolate import interp1d
import subprocess
from collections import OrderedDict
from scipy.integrate import quad
from scipy.optimize import minimize

#import json

def lc(parameters):
    
    arg = ('./disco' + ' ' + 
           str(parameters['q']) + ' ' +
           str(parameters['mu']) + ' ' +
           str(parameters['beta']) + ' ' +  
           str(parameters['u']) + ' ' + 
           str(parameters['albedo']) + ' ' +
           str(parameters['Lx']) + ' ' + 
           str(parameters['h']) + ' ' + 
           str(parameters['R']) + ' ' + 
           str(parameters['y_tilt']) + ' ' + 
           str(parameters['z_tilt']) + ' ' + 
           str(parameters['picture']) + ' ' +
           str(parameters['inclination']) + ' ' + 
           str(int(parameters['lc_num'])) + ' ' + 
           str(int(parameters['star_tiles'])) + ' ' + 
           str(int(parameters['disk_tiles'])) + ' ' +
           str(int(parameters['threads'])) + ' ' +
           str(parameters['T_disk']) + ' ' + 
           str(parameters['T_star']) + ' ' + 
           str(parameters['a']) + ' ' +
           str(parameters['y_tilt2']) + ' ' + 
           str(parameters['z_tilt2']) + ' ' +
           str(parameters['PSI_pr']) + ' ' +
           str(parameters['kappa']) + ' ' +
           str(parameters['Lx_disk']) + ' ' +
           str(parameters['spot_disk']) + ' ' +
           str(parameters['T_spot']) + ' ' +
           str(parameters['spot_beg']) + ' ' +
           str(parameters['spot_end']) + ' ' +
           str(parameters['ns_theta']) + ' ' +
           str(parameters['spot_rho_in']) + ' ' +
           str(parameters['spot_rho_out']) + ' ' +
           str(parameters['drd_phi']) + ' ' +
           str(parameters['drd_theta']) + ' ' +
           str(parameters['Lx_disk_2']) + ' ' +
           str(parameters['Lx_iso']) + ' ' +
           str(parameters['rho_in']) + ' ' +           
           str(parameters['A']) + ' ' +
           str(parameters['uniform_disk']) + ' ' +
           str(parameters['disk_flux']) + ' ' +
           str(parameters['h_warp']) + ' ' +
           str(parameters['filter']) + ' ' +
           str(parameters['lc_start']) + ' ' +
           str(parameters['lc_end'])
    )

    ##print(arg)
    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,2])
    
    x = s[:,0] + 0.5
    y = s[:,1]
    
    return x, y
    #return interp1d(x, y, kind='cubic')




def z_tilt_model(parameters):

    phase_index = float(parameters['n_data'])
    N = parameters['number_of_orbits']

    phase = (1.0/N)*phase_index

    parameters['z_tilt2'] = parameters['z_tilt'] - parameters['Z']

    parameters['PSI_pr'] = int((phase * 360.0 + parameters['deltaPsi']) % 360)

    return parameters


def residual(lmfit_parameters, parameters, data):
        
    for x in lmfit_parameters:
        parameters[str(x)] = lmfit_parameters[str(x)].value

    parameters = z_tilt_model(parameters)

    n = parameters['n_data'] 
    r = {'B':[], 'V':[]} # residuals
    
    for UBV_filter in data[n].keys():

        if data[n][UBV_filter]['flux'] == []:
            continue
        
        parameters['filter'] = UBV_filter


        x, y = lc(parameters)
        f = interp1d(x, y, kind='cubic')

        #r[UBV_filter] = data[n][UBV_filter]['flux'] - f(data[n][UBV_filter]['orbital_phase'])
        r[UBV_filter] = data[n][UBV_filter]['flux'] - f(data[n][UBV_filter]['orbit'])

    #print("Z {}\t h {}\t out {}\t in {}\t {}".format(parameters['Z'], parameters['h'], parameters['z_tilt'], parameters['z_tilt2'], (sum(np.array(r['B'])**2.0 ) + sum(np.array(r['V'])**2.0 ) )/(len(r['B']) + len(r['V'])) ))
    print("{:.1f}\t {:.1f}\t {:.1f}\t {:.1f}\t {:.1e}\t {:.5f}".format(parameters['z_tilt'], parameters['y_tilt'], parameters['y_tilt2'], parameters['Z'], parameters['Lx'], (sum(np.array(r['B'])**2.0 ) + sum(np.array(r['V'])**2.0 ) )/(len(r['B']) + len(r['V'])) ))
    
    return np.append(r['B'], r['V'])
    #return r['V']



def residual_WASP(lmfit_parameters, parameters, data):
        
    for x in lmfit_parameters:
        parameters[str(x)] = lmfit_parameters[str(x)].value

    parameters = z_tilt_model(parameters)

    n = parameters['n_data']
    
    for UBV_filter in data[n].keys():

        if data[n][UBV_filter]['flux'] == []:
            continue
        
        parameters['filter'] = UBV_filter


        x, y = lc(parameters)
        f = interp1d(x, y, kind='cubic')

        r = (data[n][UBV_filter]['flux'] - f(data[n][UBV_filter]['orbital_phase']))/np.sqrt(data[n][UBV_filter]['flux_error'])

    #print("Z {}\t h {}\t out {}\t in {}\t {}".format(parameters['Z'], parameters['h'], parameters['z_tilt'], parameters['z_tilt2'], (sum(np.array(r['B'])**2.0 ) + sum(np.array(r['V'])**2.0 ) )/(len(r['B']) + len(r['V'])) ))
    print("Z {}\t h {}\t out {}\t in {}\t {}".format(parameters['Z'], parameters['h'], parameters['z_tilt'], parameters['z_tilt2'], (sum(r**2.0)/(len(r)))))

    
    return r
    #return r['V']








pi2 = 2.0 * np.pi

dtr = (pi2/360.0)
rtd = (360.0/pi2)



def sin_epsilon(i, theta, phi, delta_phi):
    return np.cos(i) * np.cos(theta) + np.sin(i) * np.sin(theta) * np.cos(phi + delta_phi)


def fit_f(x, phi_1, phi_2, delta_phi, w):

    i = x[0]
    theta = x[1]
    
    return (np.sin(w) - sin_epsilon(i, theta, phi_1, delta_phi))**2.0 + (np.sin(-w) - sin_epsilon(i, theta, phi_2, delta_phi))**2.0



def i_and_theta_calculator(max_opening_phase, w):

    delta_phi = (-1.0) * max_opening_phase * pi2
    
    phi_1 = 0.25 * 0.05 * pi2
    phi_2 = 12.0 * 0.05 * pi2

    x0 = [pi2/4.0, pi2/20.0] # initial_guess
    
    bnds = ((80.0*dtr, 90.0*dtr), (5.0*dtr, 25.0*dtr))

    res = minimize(fit_f, x0, args=(phi_1, phi_2, delta_phi, w), method = 'TNC', bounds = bnds, tol=1e-6)

    return res.x
