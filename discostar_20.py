#from __future__ import print_function
import numpy as np
import subprocess
import datetime
from scipy.interpolate import interp1d

direct = "./LC2" # directory for LC`s
logs = "LOGS"

with open("parameters", 'r') as f:
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
1
disk_flux = p['disk_flux']


#Lx_list
#PSI_pr_list
#z_tilt_list
#y_tilt2_list
#z_tilt2_list
#T_disk_list
#y_tilt_list
#disk_flux_list
#kappa_list
#ns_theta_list


with open("./temp/phases_new_algorithm/parameters_new_data", 'r') as f:
    Lx_list, PSI_pr_list, z_tilt_list, y_tilt2_list, z_tilt2_list, T_disk_list, y_tilt_list, disk_flux_list, kappa_list, ns_theta_list = np.loadtxt(f, dtype=('float'), usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)


phase = {
    0:   '25_30.dat',
    342: '20_25.dat',
    324: '15_20.dat',
    306: '10_15.dat',
    288: '05_10.dat',
    270: '00_05.dat',
    252: '95_00.dat',
    234: '90_95.dat',
    216: '85_90.dat',
    198: '80_85.dat',
    180: '75_80.dat',
    162: '70_75.dat',
    144: '65_70.dat',
    126: '60_65.dat',
    108: '55_60.dat',
    90:  '50_55.dat',
    72:  '45_50.dat',
    54:  '40_45.dat',
    36:  '35_40.dat',
    18:  '30_35.dat',
}
 

border = 0.13
chisq_n = 1
chisq_sum = 0




for i in range(241):

    p['Lx'] = Lx_list[i] 
    p['PSI_pr'] = PSI_pr_list[i] 
    p['z_tilt'] = z_tilt_list[i] 
    p['y_tilt2'] = y_tilt2_list[i] 
    p['z_tilt2'] = z_tilt2_list[i] 
    p['T_disk'] = T_disk_list[i] 
    p['y_tilt'] = y_tilt_list[i] 
    p['disk_flux'] = disk_flux_list[i] 
    p['kappa'] = kappa_list[i] 
    p['ns_theta'] = ns_theta_list[i] 

    data_file = phase[p['z_tilt']]

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
           str(p['disk_flux'])
    )

    

    f = open('./'+direct+'/'+output_filename, "w")    
    f_logs = open('./'+logs, "a")

    subprocess.call(arg, stdout=f, shell=True)
    #subprocess.call(arg, shell=True)

    output = subprocess.check_output(arg, shell=True)
    
    s = np.array(output.split(), dtype='float').reshape([-1,2])
    
    xx = s[:,0] + 0.5
    yy = s[:,1]
    
    spline = interp1d(xx, yy, kind='cubic')

    x_short = []
    data_short = []

    for j, X in enumerate(x):
        if (X > border) or (X < 1.0 - border):
            x_short.append(x[j])
            data_short.append(data[j])

    N = len(data_short)
    chisq = np.sum(np.square(data_short - spline(x_short)))/N

    
    f_logs.write(output_filename+' '+arg+'\n')

    f.close
    f_logs.close

    print 'z_tilt', '\t', z_tilt_list[i], '\t', 'kappa', '\t', kappa_list[i], '\t', 'ns_theta', '\t', ns_theta_list[i], '\t', 'FILE', '\t', output_filename, '\t', 'chisq', '\t', chisq
        
    if chisq_n < 20:
        chisq_n = chisq_n + 1
        chisq_sum = chisq_sum + chisq
    elif chisq_n == 20:
        print 'CHISQ SUM', '\t', chisq_n, '\t', chisq_sum, '\n'
        print '\n'
        chisq_n = 1
        chisq_sum = 0
