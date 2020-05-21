#from __future__ import print_function
import numpy as np
import subprocess
import datetime

direct = "./" # directory for LC`s
logs = "LOGS"

with open("parameters", 'r') as f:
    name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    





    
p = dict(zip(name, value))   

if p['do_lc'] in ['true', 'True', 'TRUE', 'T', '1']:
    p['do_lc'] = '1'
elif p['do_lc'] in ['false', 'False', 'FALSE', 'F', '0']:
    p['do_lc'] = '0'


output_filename = "TEST.data"

    
arg = ('./disco' + ' ' + 
       p['q'] + ' ' + 
       p['mu'] + ' ' +
       p['beta'] + ' ' +
       p['u'] + ' ' +   
       p['X_albedo'] + ' ' + 
       p['T_star_polar'] + ' ' +
       p['a'] + ' ' + 
       p['inclination'] + ' ' +
       p['Lx'] + ' ' + 
       p['NS_phi'] + ' ' +
       p['NS_kappa'] + ' ' +
       p['NS_theta'] + ' ' +
       p['h_out'] + ' ' + 
       p['r_out'] + ' ' + 
       p['r_in'] + ' ' + 
       p['gamma'] + ' ' + 
       p['theta_out'] + ' ' + 
       p['phi_out'] + ' ' +
       p['theta_in'] + ' ' + 
       p['phi_in'] + ' ' +  
       p['do_lc'] + ' ' + 
       p['N_lc'] + ' ' + 
       p['N_theta'] + ' ' + 
       p['N_r'] + ' ' + 
       p['OMP_threads'] + ' ' + 
       p['filter']
)


print(arg)

    
f = open('./'+direct+'/'+output_filename, "w")    
f_logs = open('./'+logs, "a")

subprocess.call(arg, stdout=f, shell=True)
#subprocess.call(arg, shell=True)

f_logs.write(output_filename+' '+arg+'\n')

f.close
f_logs.close

print("done")

print(output_filename)
