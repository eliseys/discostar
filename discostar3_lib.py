#from __future__ import print_function
import numpy as np
from scipy.interpolate import interp1d
import subprocess
import matplotlib.pyplot as plt

def light_curve(*args):
    
    ###############################################
    # set values from parameters file
    
    with open("parameters", 'r') as f:
        parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
        
        p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}
        
    ###############################################
    ###############################################
    ###############################################
    ###############################################
    ###############################################
    # change values 

    p['picture'] = 0.0

    #print 'ARGS', args
    
    Lx = args[0]
    y_tilt2 = args[1]
    z_tilt2 = args[2]

    z_tilt = args[3]
    R = args[4]
    T_disk = args[5]

    ####
    
    p['Lx'] = Lx
    p['y_tilt2'] = y_tilt2
    p['z_tilt2'] = z_tilt2

    p['z_tilt'] = z_tilt
    p['R'] = R
    p['T_disk'] = T_disk

    ###############################################

    arg = ('./disco' + ' ' + 
           str(p['q']) + ' ' +
           str(p['mu']) + ' ' +
           str(p['beta']) + ' ' +  
           str(p['u']) + ' ' + 
           str(p['albedo']) + ' ' +
           str(p['Lx']) + ' ' + 
           str(p['h']) + ' ' + 
           str(p['R']) + ' ' + 
           str(p['y_tilt']) + ' ' + 
           str(p['z_tilt']) + ' ' + 
           str(p['picture']) + ' ' +
           str(p['inclination']) + ' ' + 
           str(int(p['lc_num'])) + ' ' + 
           str(int(p['star_tiles'])) + ' ' + 
           str(int(p['disk_tiles'])) + ' ' +
           str(int(p['threads'])) + ' ' +
           str(p['T_disk']) + ' ' + 
           str(p['T_star']) + ' ' + 
           str(p['lambda_A']) + ' ' +
           str(p['a']) + ' ' +
           str(p['y_tilt2']) + ' ' + 
           str(p['z_tilt2']) + ' ' +
           str(p['PSI_pr']) + ' ' +
           str(p['kappa']) + ' ' +
           str(p['isotrope']) + ' ' +
           str(p['Lx_disk']) + ' ' +
           str(p['spot_disk']) + ' ' +
           str(p['T_spot']) + ' ' +
           str(p['spot_beg']) + ' ' +
           str(p['spot_end']) + ' ' +
           str(p['ns_theta']) + ' ' +
           str(p['spot_rho_in']) + ' ' +
           str(p['spot_rho_out']) + ' ' +
           str(p['drd_phi']) + ' ' +
           str(p['drd_theta']) + ' ' +
           str(p['Lx_disk_2']) + ' ' +
           str(p['Lx_iso']) + ' ' +
           str(p['rho_in'])
           
    )

    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,3])
    
    x = s[:,0] + 0.5
    y = s[:,1]
    
    return interp1d(x, y, kind='cubic')


def light_curve_disk(*args):
    
    ###############################################
    # set values from parameters file
    
    with open("parameters", 'r') as f:
        parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
        
        p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}
        
    ###############################################
    ###############################################
    ###############################################
    ###############################################
    ###############################################
    # change values 

    p['picture'] = 0.0

    #print 'ARGS', args
    
    y_tilt = args[0]
    R = args[1]
    T_disk = args[2]

    #print y_tilt, y_tilt2, z_tilt
    
    p['y_tilt'] = y_tilt
    p['R'] = R
    p['T_disk'] = T_disk
    p['h'] = 0.2 * R
    
    ###############################################

    arg = ('./disco' + ' ' + 
           str(p['q']) + ' ' +
           str(p['mu']) + ' ' +
           str(p['beta']) + ' ' +  
           str(p['u']) + ' ' + 
           str(p['albedo']) + ' ' +
           str(p['Lx']) + ' ' + 
           str(p['h']) + ' ' + 
           str(p['R']) + ' ' + 
           str(p['y_tilt']) + ' ' + 
           str(p['z_tilt']) + ' ' + 
           str(p['picture']) + ' ' +
           str(p['inclination']) + ' ' + 
           str(int(p['lc_num'])) + ' ' + 
           str(int(p['star_tiles'])) + ' ' + 
           str(int(p['disk_tiles'])) + ' ' +
           str(int(p['threads'])) + ' ' +
           str(p['T_disk']) + ' ' + 
           str(p['T_star']) + ' ' + 
           str(p['lambda_A']) + ' ' +
           str(p['a']) + ' ' +
           str(p['y_tilt2']) + ' ' + 
           str(p['z_tilt2']) + ' ' +
           str(p['PSI_pr']) + ' ' +
           str(p['kappa']) + ' ' +
           str(p['isotrope']) + ' ' +
           str(p['Lx_disk']) + ' ' +
           str(p['spot_disk']) + ' ' +
           str(p['T_spot']) + ' ' +
           str(p['spot_beg']) + ' ' +
           str(p['spot_end']) + ' ' +
           str(p['ns_theta']) + ' ' +
           str(p['spot_rho_in']) + ' ' +
           str(p['spot_rho_out']) + ' ' +
           str(p['drd_phi']) + ' ' +
           str(p['drd_theta']) + ' ' +
           str(p['Lx_disk_2']) + ' ' +
           str(p['Lx_iso']) + ' ' +
           str(p['rho_in'])           
    )

    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,3])
    
    x = s[:,0] + 0.5
    y = s[:,1]
    
    return interp1d(x, y, kind='cubic')
      


def residuals(pars, result_disk, path_to_data_file):

    Lx = pars[0]
    y_tilt2 = pars[1]
    z_tilt2 = pars[2]

    y_tilt = result_disk[0]
    R = result_disk[1]
    T_disk = result_disk[2]
    
    with open(path_to_data_file, 'r') as data:
        x_obs_array, y_obs_array = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)
        
    N = len(x_obs_array)

    f = light_curve(Lx, y_tilt2, z_tilt2, y_tilt, R, T_disk)
    r = (y_obs_array - f(x_obs_array))**2
    
    return r.sum()/N



def residuals_disk(pars, border, path_to_data_file):

    y_tilt = pars[0]
    R = pars[1]
    T_disk = pars[2]

    with open(path_to_data_file, 'r') as data:
        x_obs_array, y_obs_array = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)

    x_obs_list_short = []
    y_obs_list_short = []
    
    for i, x in enumerate(x_obs_array):
        if (x < border) or (x > 1.0 - border):
            x_obs_list_short.append(x_obs_array[i])
            y_obs_list_short.append(y_obs_array[i])
        
    N = len(x_obs_list_short)

    f = light_curve_disk(y_tilt, R, T_disk)
    r = (y_obs_list_short - f(x_obs_list_short))**2

    #x = np.arange(0.0, 1.0, 0.01)
    #y = f(x)

    #plt.plot(x, y, '-', x_obs_list_short, y_obs_list_short, '.')
    #plt.show()
    
    return r.sum()/N
