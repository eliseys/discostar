import numpy as np
from scipy.interpolate import interp1d
import subprocess
from collections import OrderedDict

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
           str(parameters['lambda_A']) + ' ' +
           str(parameters['a']) + ' ' +
           str(parameters['y_tilt2']) + ' ' + 
           str(parameters['z_tilt2']) + ' ' +
           str(parameters['PSI_pr']) + ' ' +
           str(parameters['kappa']) + ' ' +
           str(parameters['isotrope']) + ' ' +
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
           str(parameters['disk_flux_B']) + ' ' +
           str(parameters['disk_flux_V']) + ' ' +
           str(parameters['h_warp']) + ' ' +
           str(parameters['UBV_filter'])
    )

    #print(arg)
    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,2])
    
    x = s[:,0] + 0.5
    y = s[:,1]
    
    return interp1d(x, y, kind='cubic')




def z_tilt_model(parameters):

    assymetry_factor = parameters['assymetry_factor'] 
    inclination = parameters['inclination']
    theta_fix = parameters['theta_fix'] 
    phase_index = float(parameters['n_data']) 
    
    delta_z_tilt_dot = (180.0/np.pi) * 0.05 * 2.0 * np.pi * (2.0 * np.pi * assymetry_factor - np.arccos(-np.tan((90.0 - inclination) * np.pi / 180.0)/np.tan(theta_fix * np.pi/180.0)))/np.sin(2.0*np.pi*assymetry_factor)

    #xdelta_z_tilt_dot = 8.0
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #if delta_z_tilt_dot < 0.0:
    #    continue
    
    phase = (1.0/20.0)*phase_index

    parameters['z_tilt'] = ((phase - (delta_z_tilt_dot/18.0) * np.sin((phase - 0.1)*2*np.pi)/(2*np.pi)) - 0.1) * 360.0 % 360.0
    
    parameters['z_tilt2'] = parameters['z_tilt'] + parameters['Z']

    parameters['PSI_pr'] = (phase * 360.0 + parameters['deltaPsi']) % 360

    return parameters




def residual(lmfit_parameters, parameters, data):
        
    for x in lmfit_parameters:
        parameters[str(x)] = lmfit_parameters[str(x)].value

    parameters = z_tilt_model(parameters)
        
    n = parameters['n_data']

    r = OrderedDict()
    
    for UBV_filter in data[n].keys():

        if data[n][UBV_filter]['flux'] == []:
            continue
        
        parameters['UBV_filter'] = UBV_filter
        f = lc(parameters)
        r[UBV_filter] = data[n][UBV_filter]['flux'] - f(data[n][UBV_filter]['orbit'])

    #print('{:8f}'.format(np.sum(np.square(data_short - f(x_short)))/N), '\t', '{:8f}'.format(y_tilt),  '\t','{:8f}'.format(y_tilt2), '\t', '{:E}'.format(Lx), '\t', '{:E}'.format(disk_flux))

    print("Res = {}".format(sum(r['B']**2.0)/len(r['B'])))
    
    return r['B']

