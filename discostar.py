#from __future__ import print_function
import numpy as np
import subprocess
import datetime

direct = "./LC2" # directory for LC`s
logs = "LOGS"

with open("parameters", 'r') as f:
    parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    



p = {}

for i in range(len(parameter_name)):
    if parameter_name[i] != 'UBV_filter':
        p[parameter_name[i]] = float(value[i])
    else:
        p[parameter_name[i]] = str(value[i])


    


#p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}




z_tilt = p['z_tilt']
Lx = p['Lx']
y_tilt = p['y_tilt']
y_tilt2 = p['y_tilt2']
z_tilt2 = p['z_tilt2']
T_disk = p['T_disk']
h = p['h']
R = p['R']

picture = p['picture']

inclination = p['inclination']
q = p['q']
mu = p['mu']
beta = p['beta']
u = p['u']
albedo = p['albedo']
lc_num = p['lc_num']
star_tiles = p['star_tiles']
disk_tiles = p['disk_tiles']
threads = p['threads']
T_star = p['T_star']
lambda_A = p['lambda_A']
a = p['a']

PSI_pr = p['PSI_pr']
kappa = p['kappa'];

isotrope = p['isotrope']

Lx_disk = p['Lx_disk']

spot_disk = p['spot_disk']

T_spot = p['T_spot']
spot_beg = p['spot_beg']
spot_end = p['spot_end']

ns_theta = p['ns_theta']

spot_rho_in = p['spot_rho_in']
spot_rho_out = p['spot_rho_out']

drd_phi = p['drd_phi']
drd_theta = p['drd_theta']

Lx_disk_2 = p['Lx_disk_2']
Lx_iso = p['Lx_iso']

rho_in = p['rho_in']
A = p['A']

uniform_disk = p['uniform_disk']

disk_flux = p['disk_flux']

h_warp = p['h_warp']


if picture == 0:
    output_filename = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")+".data"        
elif picture == 1:
    output_filename = 'VIEW.data'
    
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
       str(p['h_warp']) + ' ' +
       str(p['UBV_filter'])

)


if picture == 0:
    print 'discostar calculates the lightcurve ...'
elif picture == 1:
    print 'discostar draws the picture ...'
    

f = open('./'+direct+'/'+output_filename, "w")    
f_logs = open('./'+logs, "a")

subprocess.call(arg, stdout=f, shell=True)
#subprocess.call(arg, shell=True)

f_logs.write(output_filename+' '+arg+'\n')

f.close
f_logs.close

#print 'done'

if picture == 0:
    print 'lightcurve written to', output_filename
        
elif picture == 1:
    print 'data for picture written to', output_filename
