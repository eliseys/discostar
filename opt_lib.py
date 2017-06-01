import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
import subprocess

with open("parameters", 'r') as f:
    parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    
p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}

# default parameters
z_tilt_initial = p['z_tilt']
z_tilt = p['z_tilt']  
Lx = p['Lx']
h = p['h']
R = p['R']
y_tilt = p['y_tilt']
z_tilt = p['z_tilt']
b = p['b']
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
T_disk = p['T_disk']
T_star = p['T_star']
lambda_A = p['lambda_A']
a = p['a']
y_tilt2 = p['y_tilt2']
z_tilt2 = p['z_tilt2'] 

# changed parameters for optimisation procedure
lc_num = 20

direct = "LC" # directory for light curves

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

def flx(Lx,
        T_disk,
        y_tilt,
        y_tilt2,
        z_tilt2,
        z_tilt):
    
# the only difference with flx_xy is the argument order
    
    arg = ['./disco', 
           str(q),
           str(mu),
           str(beta) ,  
           str(u) , 
           str(albedo) ,
           str(Lx) , 
           str(h) , 
           str(R) , 
           str(y_tilt) , 
           str(z_tilt) , 
           str(b) ,
           str(inclination) , 
           str(lc_num) , 
           str(star_tiles) , 
           str(disk_tiles) ,
           str(threads) ,
           str(T_disk) , 
           str(T_star) , 
           str(lambda_A) ,
           str(a) ,
           str(y_tilt2) , 
           str(z_tilt2)
    ]

    output = subprocess.check_output(arg)

    s = np.array(output.split(), dtype='float').reshape([-1,2])
    
    x = s[:,0] + 0.5
    y = s[:,1]
    
    return interp1d(x, y, kind='cubic')


y_obs = []
x_obs = []

inguess_LT = np.array([Lx, T_disk])
inguess_yz = np.array([y_tilt, y_tilt2, z_tilt2])

def fit_LT(normalized_parameters, y_tilt, y_tilt2, z_tilt2, z_tilt):
    
    #
    # the only difference with fit_yz is the argument order
    # fit_LT(Lx, T_disk, y_tilt, y_tilt2, z_tilt2, z_tilt)
    #
    
    parameters = inguess_LT * normalized_parameters

    Lx = parameters[0]
    T_disk = parameters[1]
    
    f = flx(Lx, T_disk, y_tilt, y_tilt2, z_tilt2, z_tilt)
    r_2 = (y_obs - f(x_obs))**2
    return r_2.sum()


def fit_yz(normalized_parameters, Lx, T_disk, z_tilt):

    #
    # the only difference with fit_LT is the argument order
    # fit_yz(y_tilt, y_tilt2, z_tilt2, Lx, T_disk, z_tilt)
    #

    parameters = inguess_yz * normalized_parameters

    y_tilt = parameters[0]
    y_tilt2 = parameters[1]
    z_tilt2 = parameters[2]
    
    f = flx(Lx, T_disk, y_tilt, y_tilt2, z_tilt2, z_tilt)
    r_2 = (y_obs - f(x_obs))**2
    return r_2.sum()
