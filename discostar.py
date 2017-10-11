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


if picture == 0:
    #output = 'LC.data'
    output = 'LC_Lx_{Lx}_{Lx_disk}_NS_{PSI_pr}_{kappa}_{ns_theta}_h_R_{h}_{R}_tilt_{y_tilt}_{z_tilt}_{y_tilt2}_{z_tilt2}_Td_{T_disk}_spot_{spot_disk}_{T_spot}_{spot_beg}_{spot_end}_drd_{drd_phi}_{drd_theta}_i_{inclination}.data'
    if (isotrope == 0 and spot_disk != 0):
        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, PSI_pr=PSI_pr, kappa=kappa, ns_theta=ns_theta, h=h, R=R, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot=T_spot, spot_beg=spot_beg, spot_end=spot_end, drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
    elif (isotrope == 1 and spot_disk != 0):
        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, PSI_pr='', kappa='', ns_theta='', h=h, R=R, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot=T_spot, spot_beg=spot_beg, spot_end=spot_end, drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
    elif (isotrope == 1 and spot_disk == 0):
        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, PSI_pr='', kappa='', ns_theta='', h=h, R=R, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot='', spot_beg='', spot_end='', drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
    elif (isotrope == 0 and spot_disk == 0):
        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, PSI_pr=PSI_pr, kappa=kappa, ns_theta=ns_theta, h=h, R=R, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot='', spot_beg='', spot_end='', drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
        
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
       str(p['drd_theta'])
)


if picture == 0:
    print 'discostar calculates the lightcurve ...'
elif picture == 1:
    print 'discostar draws the picture ...'
    
#print '-----------------------------------------------------------------------------------------'
#print '|', 'q', '\t', p['q'], '\t', '|', 'albedo', '\t', p['albedo'],'\t', '|', '     ', '\t', '\t', '|', 'lc_num', '\t', p['lc_num'],'\t', '|'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'mu', '\t', p['mu'],'\t', '|', 'Lx','\t', '\t', p['Lx'], '\t', '|', '     ', '\t', '\t', '|', 'star_tiles', '\t', p['star_tiles'],'\t', '|'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'beta', '\t', p['beta'],'\t', '|', 'h','\t', '\t', p['h'], '\t', '|', ' ','\t', '\t', '\t', '|', 'disk_tiles', '\t', p['disk_tiles'],'\t', '|'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'u', '\t', p['u'],'\t', '|', 'R','\t', '\t', p['R'], '\t', '|', '     ', '\t', '\t', '|', 'threads', '\t', p['threads'],'\t', '|'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'T_star', '\t', p['T_star'], '\t', '|', 'T_disk', '\t', p['T_disk'], '\t', '|', 'lambda_A', '\t', p['lambda_A'], '\t', '\t', '|'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'PSI_pr', '\t', p['PSI_pr'], '\t',  '|', 'kappa', '\t', p['kappa'], '\t'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'y_tilt', '\t', p['y_tilt'],'\t', '|', 'y_tilt2', '\t', p['y_tilt2'],'\t'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'z_tilt', '\t', p['z_tilt'],'\t', '|', 'z_tilt2', '\t', p['z_tilt2'],'\t'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'inclination', '\t', p['inclination'],'\t'
#print '-----------------------------------------------------------------------------------------'
#print '|', 'Lx_disk', '\t', p['Lx_disk'],'\t'
#print '-----------------------------------------------------------------------------------------'


f = open('./'+direct+'/'+output_filename, "w")    
subprocess.call(arg, stdout=f, shell=True)
#subprocess.call(arg, shell=True)

f.close

print 'done'

if picture == 0:
    print 'lightcurve written to', output_filename
    print '\n'
    print '\n'
    if (spot_disk != 0 and isotrope == 0):
        print '"../LC2/%s" @legend "%3.0f (%3.0f %2.0f %2.0f) (%3.0f %2.0f %2.0f) hR=%1.2f (%1.1e %1.1e) Td=%5.0f SPOT %d",\\' % (output_filename,z_tilt,PSI_pr,kappa,ns_theta,z_tilt2,y_tilt,y_tilt2,h/R,Lx,Lx_disk,T_disk,spot_disk)
    elif (spot_disk != 0 and isotrope == 1):
        print '"../LC2/%s" @legend "%3.0f ISO (%3.0f %2.0f %2.0f) hR=%1.2f (%1.1e %1.1e) Td=%5.0f SPOT %d",\\' % (output_filename,z_tilt,z_tilt2,y_tilt,y_tilt2,h/R,Lx,Lx_disk,T_disk,spot_disk)
    elif (spot_disk == 0 and isotrope == 1):
        print '"../LC2/%s" @legend "%3.0f ISO (%3.0f %2.0f %2.0f) hR=%1.2f (%1.1e %1.1e) Td=%5.0f NO SPOT",\\' % (output_filename,z_tilt,z_tilt2,y_tilt,y_tilt2,h/R,Lx,Lx_disk,T_disk)
    elif (spot_disk == 0 and isotrope == 0):
        print '"../LC2/%s" @legend "%3.0f (%3.0f %2.0f %2.0f) (%3.0f %2.0f %2.0f) hR=%1.2f (%1.1e %1.1e) Td=%5.0f NO SPOT",\\' % (output_filename,z_tilt,PSI_pr,kappa,ns_theta,z_tilt2,y_tilt,y_tilt2,h/R,Lx,Lx_disk,T_disk)   
    print '\n'
    print '\n'

        
elif picture == 1:
    print 'data for picture written to', output_filename
    





