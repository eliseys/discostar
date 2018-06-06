#from __future__ import print_function
import numpy as np
import subprocess
import datetime
from scipy.interpolate import interp1d

direct = "./LC2" # directory for LC`s
logs = "LOGS"

with open("parameters_copy", 'r') as f:
    parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    
p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}

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

with open("./temp/phases_new_algorithm/output_i87_z_tilt2_const_set_kappa_plus15_minus15_0", 'r') as f:
    #phase_list, rchisq_list, Lx_list, PSI_pr_list, z_tilt_list, y_tilt_list, z_tilt2_list, y_tilt2_list, kappa_list, ns_theta_list  = np.loadtxt(f, dtype=('float'), usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
    #phase_list, rchisq_list, Lx_list, PSI_pr_list, z_tilt_list, y_tilt_list, y_tilt2_list, z_tilt2_list, T_disk_list, kappa_list, ns_theta_list  = np.loadtxt(f, dtype=('float'), usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
    phase_list, rchisq_list, Lx_list, PSI_pr_list, z_tilt_list, y_tilt_list, z_tilt2_list, T_disk_list, y_tilt2_list, kappa_list, ns_theta_list = np.loadtxt(f, dtype=('float'), usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)


phase = {
    0.00: '00_05.dat',
    0.05: '05_10.dat',
    0.10: '10_15.dat',
    0.15: '15_20.dat',
    0.20: '20_25.dat',
    0.25: '25_30.dat',
    0.30: '30_35.dat',
    0.35: '35_40.dat',
    0.40: '40_45.dat',
    0.45: '45_50.dat',
    0.50: '50_55.dat',
    0.55: '55_60.dat',
    0.60: '60_65.dat',
    0.65: '65_70.dat',
    0.70: '70_75.dat',
    0.75: '75_80.dat',
    0.80: '80_85.dat',
    0.85: '85_90.dat',
    0.90: '90_95.dat',
    0.95: '95_00.dat',  
}


border = 0.13
chisq_n = 1
chisq_sum = 0


for i in range(302):


    phi = phase_list[i]
    rchisq = rchisq_list[i]
    
    p['Lx'] = Lx_list[i] 
    p['PSI_pr'] = PSI_pr_list[i]
   
    p['z_tilt'] = z_tilt_list[i] 

    p['y_tilt2'] = y_tilt2_list[i] 
    p['z_tilt2'] = z_tilt2_list[i] 
    #p['T_disk'] = T_disk_list[i] 
    p['y_tilt'] = y_tilt_list[i] 
    p['kappa'] = kappa_list[i] 
    p['ns_theta'] = ns_theta_list[i] 

    data_file = phase[phi]

    path_to_data_file = './OBSERVED_DATA/'+data_file

    with open(path_to_data_file, 'r') as data:
        x, data = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)
    
    output_filename = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")+".data"
    
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

    

    f = open('./'+direct+'/'+output_filename, "w")    
    f_logs = open('./'+logs, "a")

    subprocess.call(arg, stdout=f, shell=True)
    #subprocess.call(arg, shell=True)

    output = subprocess.check_output(arg, shell=True)
    
    #s = np.array(output.split(), dtype='float').reshape([-1,2])
    
    #xx = s[:,0] + 0.5
    #yy = s[:,1]
    
    #spline = interp1d(xx, yy, kind='cubic')

    #x_short = []
    #data_short = []

    #for j, X in enumerate(x):
    #    if (X > border) or (X < 1.0 - border):
    #        x_short.append(x[j])
    #        data_short.append(data[j])

    #N = len(data_short)
    #chisq = np.sum(np.square(data_short - spline(x_short)))/N

    
    f_logs.write(output_filename+' '+arg+'\n')

    f.close
    f_logs.close

    print 'phi', '\t', phi, '\t', '\t', 'FILE', '\t', output_filename, '\t', 'chisq', '\t',  '{:6.5f}'.format(rchisq)
        
#    if chisq_n < 20:
#        chisq_n = chisq_n + 1
#        chisq_sum = chisq_sum + chisq
#    elif chisq_n == 20:
#        print 'CHISQ SUM', '\t', chisq_n, '\t', chisq_sum, '\n'
#        print '\n'
#        chisq_n = 1
#        chisq_sum = 0
