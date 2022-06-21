#from __future__ import print_function
import numpy as np
import subprocess
import datetime

direct = "./" # directory for LC`s
logs = "LOGS"

with open("parameters", 'r') as f:
    name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    
    
p = dict(zip(name, value))   



def boolifier(p, p_name):
    try:
        if p[p_name] in ['true', 'True', 'TRUE', 'T', '1']:
            p[p_name] = '1'
        elif p[p_name] in ['false', 'False', 'FALSE', 'F', '0']:
            p[p_name] = '0'
        else:
            raise ValueError(p[p_name])
    except ValueError:
        raise

boolifier(p, 'do_lc')
boolifier(p, 'do_corona')


# filter parameters check
# try:
#     if (p['filter'] in ['B', 'V', 'WASP', 'vis']):
#         pass
#     else:
#         float(p['filter'])
# except ValueError:
#     print("filter must be either float number or one of the following: 'B', 'V', 'WASP', 'vis'")
#     raise



# convert degrees to radians
p['inclination'] = str(float(p['inclination']) * (np.pi/180.0))
p['NS_phi'] = str(float(p['NS_phi']) * (np.pi/180.0))
p['NS_kappa'] = str(float(p['NS_kappa']) * (np.pi/180.0))
p['NS_theta'] = str(float(p['NS_theta']) * (np.pi/180.0))
p['theta_out'] = str(float(p['theta_out']) * (np.pi/180.0))
p['phi_out'] = str(float(p['phi_out']) * (np.pi/180.0))
p['theta_in']  = str(float(p['theta_in']) * (np.pi/180.0))
p['phi_in']  = str(float(p['phi_in']) * (np.pi/180.0))


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
       p['do_corona'] + ' ' + 
       p['N_corona'] + ' ' + 
       p['h_corona'] + ' ' + 
       p['rs_corona']  + ' ' +
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
