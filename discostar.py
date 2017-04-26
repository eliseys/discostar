import numpy as np
import math
import subprocess

direct = "LC" # directory for LC`s

with open("parameters", 'r') as f:
    parameter_name, value = np.loadtxt(f, dtype=('str'), usecols=(0,1), unpack=True)
    
p = {parameter_name[i]: float(value[i]) for i in range(len(parameter_name))}

# parameters we are searching for
range_Lx = [20.0, 180.0]
range_h = [0.001, 0.1]
range_R = [0.05, 0.4]
range_y_tilt = [0.0, 50.0]
range_b = [1.0, 50.0]
range_inclination = [70.0, 90.0]

# number of nodes of the array
K = 1

for i in range(K+1):

    Lx = range_Lx[0] + i * (range_Lx[-1] - range_Lx[0])/K 
    h = range_h[0] + i * (range_h[-1] - range_h[0])/K 
    R = range_R[0] + i * (range_R[-1] - range_R[0])/K 
    y_tilt = range_y_tilt[0] + i * (range_y_tilt[-1] - range_y_tilt[0])/K 
    b = range_b[0] + i * (range_b[-1] - range_b[0])/K 
    inclination = range_inclination[0] + i * (range_inclination[-1] - range_inclination[0])/K

    p['Lx'] = Lx
    p['h'] = h
    p['R'] = R
    p['y_tilt'] = y_tilt
    p['b'] = b
    p['inclination'] = inclination
    
    q = p['q']
    mu = p['mu']
    beta = p['beta']
    u = p['u']
    albedo = p['albedo']
    lc_num = p['lc_num']
    star_tiles = p['star_tiles']
    disk_tiles = p['disk_tiles']
    threads = p['threads']
    z_tilt = p['z_tilt']

    output = 'LC_{q}_{mu}_{beta}_{u}_{albedo}_{Lx}_{h}_{R}_{y_tilt}_{z_tilt}_{b}_{inclination}_{lc_num}_{star_tiles}_{disk_tiles}_{threads}.data'
    output_filename = output.format(q=q, mu=mu, beta=beta, u=u, albedo=albedo, Lx=Lx, h=h, R=R, y_tilt=y_tilt, z_tilt=z_tilt, b=b, inclination=inclination, lc_num=lc_num, star_tiles=star_tiles, disk_tiles=disk_tiles, threads=threads)
 
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
           str(int(p['threads']))
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

    f = open('./'+direct+'/'+output_filename, "w")    
        
    subprocess.call(arg, stdout=f, shell=True)
        
    f.close
