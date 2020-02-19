#from __future__ import print_function
import numpy as np
import subprocess
import datetime

direct = "./" # directory for LC`s
logs = "LOGS"

with open("parameters", 'r') as f:
    parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    

p = {}

for i in range(len(parameter_name)):
    if parameter_name[i] != 'filter': # this is because UBV_filter is the only string parameter
        p[parameter_name[i]] = float(value[i])
    else:
        p[parameter_name[i]] = str(value[i])


if p['picture'] == 0:
    output_filename = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")+".data"        
elif p['picture'] == 1:
    output_filename = 'VIEW.data'
    
arg = ('./disco' + ' ' + 
       str(p['q']) + ' ' +
       str(p['mu']) + ' ' +
       str(p['beta']) + ' ' +  
       str(p['u']) + ' ' + 
       str(p['albedo']) + ' ' +
       str(p['Lx_noniso']) + ' ' + 
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
       str(p['a']) + ' ' +
       str(p['y_tilt2']) + ' ' + 
       str(p['z_tilt2']) + ' ' +
       str(p['PSI_pr']) + ' ' +
       str(p['kappa']) + ' ' +
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
       str(p['h_warp']) + ' ' +
       str(p['filter']) + ' ' +
       str(p['lc_start']) + ' ' +
       str(p['lc_end'])
)


print(arg)


if p['picture'] == 0:
    print 'discostar calculates the lightcurve ...'
elif p['picture'] == 1:
    print 'discostar draws the picture ...'
    

f = open('./'+direct+'/'+output_filename, "w")    
f_logs = open('./'+logs, "a")

subprocess.call(arg, stdout=f, shell=True)
#subprocess.call(arg, shell=True)

f_logs.write(output_filename+' '+arg+'\n')

f.close
f_logs.close

#print 'done'

if p['picture'] == 0:
    print 'lightcurve written to', output_filename
        
elif p['picture'] == 1:
    print 'data for picture written to', output_filename
