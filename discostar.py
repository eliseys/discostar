#from __future__ import print_function
import numpy as np
import subprocess

direct = "./LC2" # directory for LC`s

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

#output = 'LC_{q}_{mu}_{beta}_{u}_{albedo}_{Lx}_{h}_{R}_{y_tilt}_{z_tilt}_{b}_{inclination}_{lc_num}_{star_tiles}_{disk_tiles}_{threads}_{T_disk}_{T_star}_{lambda_A}_{y_tilt2}_{z_tilt2}.data'
#output_filename = output.format(q=q, mu=mu, beta=beta, u=u, albedo=albedo, Lx=Lx, h=h, R=R, y_tilt=y_tilt, z_tilt=z_tilt, b=b, inclination=inclination, lc_num=lc_num, star_tiles=star_tiles, disk_tiles=disk_tiles, threads=threads, T_disk = T_disk, T_star = T_star, lambda_A = lambda_A, y_tilt2 = y_tilt2, z_tilt2 = z_tilt2)

if picture == 0:
    output = 'LC_{z_tilt}.data'
elif picture == 1:
    output = 'LOOK_{z_tilt}.data'

output_filename = output.format(z_tilt = z_tilt)

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
       str(p['isotrope'])
)


print 'discostar calculates the light curve ... '
print '-----------------------------------------------------------------------------------------'
print '|', 'q', '\t', p['q'], '\t', '|', 'albedo', '\t', p['albedo'],'\t', '|', '     ', '\t', '\t', '|', 'lc_num', '\t', p['lc_num'],'\t', '|'
print '-----------------------------------------------------------------------------------------'
print '|', 'mu', '\t', p['mu'],'\t', '|', 'Lx','\t', '\t', p['Lx'], '\t', '|', '     ', '\t', '\t', '|', 'star_tiles', '\t', p['star_tiles'],'\t', '|'
print '-----------------------------------------------------------------------------------------'
print '|', 'beta', '\t', p['beta'],'\t', '|', 'h','\t', '\t', p['h'], '\t', '|', ' ','\t', '\t', '\t', '|', 'disk_tiles', '\t', p['disk_tiles'],'\t', '|'
print '-----------------------------------------------------------------------------------------'
print '|', 'u', '\t', p['u'],'\t', '|', 'R','\t', '\t', p['R'], '\t', '|', '     ', '\t', '\t', '|', 'threads', '\t', p['threads'],'\t', '|'
print '-----------------------------------------------------------------------------------------'
print '|', 'T_star', '\t', p['T_star'], '\t', '|', 'T_disk', '\t', p['T_disk'], '\t', '|', 'lambda_A', '\t', p['lambda_A'], '\t', '\t', '|'
print '-----------------------------------------------------------------------------------------'
print '|', 'PSI_pr', '\t', p['PSI_pr'], '\t',  '|', 'kappa', '\t', p['kappa'], '\t'
print '-----------------------------------------------------------------------------------------'
print '|', 'y_tilt', '\t', p['y_tilt'],'\t', '|', 'y_tilt2', '\t', p['y_tilt2'],'\t'
print '-----------------------------------------------------------------------------------------'
print '|', 'z_tilt', '\t', p['z_tilt'],'\t', '|', 'z_tilt2', '\t', p['z_tilt2'],'\t'
print '-----------------------------------------------------------------------------------------'
print '|', 'inclination', '\t', p['inclination'],'\t'
print '-----------------------------------------------------------------------------------------'


f = open('./'+direct+'/'+output_filename, "w")    
subprocess.call(arg, stdout=f, shell=True)
#subprocess.call(arg, shell=True)

f.close


