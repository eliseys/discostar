import numpy as np
from scipy.interpolate import interp1d
import subprocess

def lc(*args):
    
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

    y_tilt = args[3]
    #R = args[4]
    T_disk = args[4]

    z_tilt = args[5]
    PSI_pr = args[6]
    
    ####

    p['z_tilt'] = z_tilt
    p['PSI_pr'] = PSI_pr
    
    p['Lx'] = Lx
    p['y_tilt2'] = y_tilt2
    p['z_tilt2'] = z_tilt2

    p['y_tilt'] = y_tilt
    #p['R'] = R
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
           str(p['rho_in']) + ' ' +
           str(p['A']) + ' ' +
           str(p['uniform_disk']) + ' ' +
           str(p['h_warp'])
    )

    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,2])
    
    x = s[:,0] + 0.5
    y = s[:,1]
    
    return interp1d(x, y, kind='cubic')


def lc_disk(*args):
    
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

    z_tilt = args[0]
    PSI_pr = args[1]

    Lx = args[2]
    y_tilt = args[3]
    y_tilt2 = args[4]
    z_tilt2 = args[5]
    disk_flux = args[6]
    T_disk = args[7]
    
    ns_theta = args[8]
    kappa = args[9]
    
    ###


    p['Lx'] = Lx

    p['z_tilt'] = z_tilt
    p['z_tilt2'] = z_tilt2
    
    p['PSI_pr'] = PSI_pr
    
    p['y_tilt'] = y_tilt
    p['y_tilt2'] = y_tilt2

    p['disk_flux'] = disk_flux
    
    p['T_disk'] = T_disk
    
    p['ns_theta'] = ns_theta
    p['kappa'] = kappa


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
           str(p['rho_in']) + ' ' +           
           str(p['A']) + ' ' +
           str(p['uniform_disk']) + ' ' +
           str(p['disk_flux']) + ' ' +
           str(p['h_warp'])
    )

    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,6])
    
    x = s[:,0] + 0.5
    y = s[:,1]

    return interp1d(x, y, kind='cubic')

