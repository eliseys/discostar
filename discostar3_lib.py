#from __future__ import print_function
import numpy as np
from scipy.interpolate import interp1d
import subprocess

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
    
    y_tilt = args[0]
    y_tilt2 = args[1]
    z_tilt = args[2]

    #print y_tilt, y_tilt2, z_tilt
    
    p['y_tilt'] = y_tilt
    p['y_tilt2'] = y_tilt2
    p['z_tilt'] = z_tilt

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
           str(p['Lx_iso'])
           
    )

    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,3])
    
    x = s[:,0] + 0.5
    y = s[:,1]
    
    return interp1d(x, y, kind='cubic')
        
