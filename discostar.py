import numpy as np
import math
import subprocess

direct = "LC" # directory for LC`s

with open("parameters", 'r') as f:
    parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    
p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}

# parameters we are searching for
#range_Lx = [100.0]
#range_h = [0.023]
#range_R = [0.266]
#range_y_tilt = [30.0]
range_z_tilt = [0.0, 180.0]
#range_b = [1.0]
#range_inclination = [90.0]

# number of nodes of the array
K = 1

for i1 in range(K+1):

    z_tilt = range_z_tilt[0] + i1 * (range_z_tilt[-1] - range_z_tilt[0])/K 

    p['z_tilt'] = z_tilt
    
    Lx = p['Lx']
    h = p['h']
    R = p['R']
    y_tilt = p['y_tilt']
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
    
    output = 'LC_{q}_{mu}_{beta}_{u}_{albedo}_{Lx}_{h}_{R}_{y_tilt}_{z_tilt}_{b}_{inclination}_{lc_num}_{star_tiles}_{disk_tiles}_{threads}_{T_disk}_{T_star}_{lambda_A}.data'
    output_filename = output.format(q=q, mu=mu, beta=beta, u=u, albedo=albedo, Lx=Lx, h=h, R=R, y_tilt=y_tilt, z_tilt=z_tilt, b=b, inclination=inclination, lc_num=lc_num, star_tiles=star_tiles, disk_tiles=disk_tiles, threads=threads, T_disk = T_disk, T_star = T_star, lambda_A = lambda_A)
 
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
           str(p['b']) + ' ' +
           str(p['inclination']) + ' ' + 
           str(int(p['lc_num'])) + ' ' + 
           str(int(p['star_tiles'])) + ' ' + 
           str(int(p['disk_tiles'])) + ' ' +
           str(int(p['threads'])) + ' ' +
           str(p['T_disk']) + ' ' + 
           str(p['T_star']) + ' ' + 
           str(p['lambda_A']) + ' ' +
           str(p['a'])
    )

                        
    print 'discostar calculates the light curve ... '
    print '-----------------------------------------------------------------------------------------'
    print '|', 'q', '\t', p['q'], '\t', '|', 'albedo', '\t', p['albedo'],'\t', '|', 'y_tilt', '\t', p['y_tilt'],'\t', '|', 'lc_num', '\t', p['lc_num'],'\t', '|'
    print '-----------------------------------------------------------------------------------------'
    print '|', 'mu', '\t', p['mu'],'\t', '|', 'Lx','\t', '\t', p['Lx'], '\t', '|', 'z_tilt', '\t', p['z_tilt'],'\t', '|', 'star_tiles', '\t', p['star_tiles'],'\t', '|'
    print '-----------------------------------------------------------------------------------------'
    print '|', 'beta', '\t', p['beta'],'\t', '|', 'h','\t', '\t', p['h'], '\t', '|', 'b','\t', '\t', p['b'],'\t', '|', 'disk_tiles', '\t', p['disk_tiles'],'\t', '|'
    print '-----------------------------------------------------------------------------------------'
    print '|', 'u', '\t', p['u'],'\t', '|', 'R','\t', '\t', p['R'], '\t', '|', 'inclination', '\t', p['inclination'],'\t', '|', 'threads', '\t', p['threads'],'\t', '|'
    print '-----------------------------------------------------------------------------------------'
    print '|', 'T_star', '\t', p['T_star'], '\t', '|', 'T_disk', '\t', p['T_disk'], '\t', '|', 'lambda_A', '\t', p['lambda_A'], '\t', '\t', '\t', '|'
    print '-----------------------------------------------------------------------------------------'
    print '|', 'a', '\t', '\t', p['a'], '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '|'
    print '-----------------------------------------------------------------------------------------'



    
    f = open('./'+direct+'/'+output_filename, "w")    
                        
    subprocess.call(arg, stdout=f, shell=True)
    #subprocess.call(arg, shell=True)
       
    f.close
                        
