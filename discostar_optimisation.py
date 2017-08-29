#from __future__ import print_function

from scipy.interpolate import interp1d
import numpy as np
import subprocess

directory = "./LC2" # directory for LC`s
data_directory = './OBSERVED_DATA/'
data_files = {
    36: '00_05.dat',	
    18: '05_10.dat',  
    0: '10_15.dat',
    342: '15_20.dat', 
    324: '20_25.dat',
    306: '25_30.dat',
    288: '30_35.dat',
    270: '35_40.dat',
    252: '40_45.dat',
    234: '45_50.dat',
    216: '50_55.dat',
    198: '55_60.dat',
    180: '60_65.dat',
    162: '65_70.dat',
    144: '70_75.dat',
    126: '75_80.dat',
    108: '80_85.dat',
    90: '85_90.dat',
    72: '90_95.dat',
    54: '95_00.dat'
}

with open("parameters", 'r') as f:
    parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    
p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}

###########################################################
###########################################################
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
###########################################################
###########################################################

#output_filename = 'discostar_optimisation_logs'

def flx(kappa):    
    arg = ['./disco',
           str(p['q']),
           str(p['mu']),
           str(p['beta']), 
           str(p['u']),
           str(p['albedo']),
           str(p['Lx']),
           str(p['h']),
           str(p['R']),
           str(p['y_tilt']),
           str(p['z_tilt']),
           str(p['picture']),
           str(p['inclination']),
           str(int(p['lc_num'])),
           str(int(p['star_tiles'])),
           str(int(p['disk_tiles'])),
           str(int(p['threads'])),
           str(p['T_disk']),
           str(p['T_star']),
           str(p['lambda_A']),
           str(p['a']),
           str(p['y_tilt2']),
           str(p['z_tilt2']),
           str(p['PSI_pr']),
           str(p['kappa']),
           str(p['isotrope'])
    ]
     
    output = subprocess.check_output(arg)
    s = np.array(output.split(), dtype='float').reshape([-1,2])
    x = s[:,0] + 0.5
    y = s[:,1]
    return interp1d(x, y, kind='cubic')



#with open(data_directory+data_files[z_tilt], 'r') as data:
#    opt_lib.x_obs, opt_lib.y_obs = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)


result = flx(30.0)
