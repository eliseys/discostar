#from __future__ import print_function
import numpy as np
import subprocess
import datetime



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

disk_flux = p['disk_flux']



if picture == 0:

    output_filename = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")+".data"

    #output = 'LC.data'
#    output = 'LC_Lx_{Lx}_{Lx_disk}_{Lx_disk_2}_{Lx_iso}_NS_{PSI_pr}_{kappa}_{ns_theta}_h_{h}_R_{R}_rho_in_{rho_in}_A_{A}_ud_{uniform_disk}_tilt_{y_tilt}_{z_tilt}_{y_tilt2}_{z_tilt2}_Td_{T_disk}_spot_{spot_disk}_{T_spot}_{spot_beg}_{spot_end}_drd_{drd_phi}_{drd_theta}_i_{inclination}.data'
#    if (isotrope == 0 and spot_disk != 0):
#        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, Lx_disk_2=Lx_disk_2, Lx_iso=Lx_iso, PSI_pr=PSI_pr, kappa=kappa, ns_theta=ns_theta, h=h, R=R, rho_in=rho_in, A=A, uniform_disk=uniform_disk, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot=T_spot, spot_beg=spot_beg, spot_end=spot_end, drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
#    elif (isotrope == 1 and spot_disk != 0):
#        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, Lx_disk_2=Lx_disk_2, Lx_iso=Lx_iso, PSI_pr='', kappa='', ns_theta='', h=h, R=R, rho_in=rho_in,  A=A, uniform_disk=uniform_disk, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot=T_spot, spot_beg=spot_beg, spot_end=spot_end, drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
#    elif (isotrope == 1 and spot_disk == 0):
#        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, Lx_disk_2=Lx_disk_2, Lx_iso=Lx_iso, PSI_pr='', kappa='', ns_theta='', h=h, R=R, rho_in=rho_in, A=A, uniform_disk=uniform_disk, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot='', spot_beg='', spot_end='', drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
#    elif (isotrope == 0 and spot_disk == 0):
#        output_filename = output.format(Lx=Lx, Lx_disk=Lx_disk, Lx_disk_2=Lx_disk_2, Lx_iso=Lx_iso, PSI_pr=PSI_pr, kappa=kappa, ns_theta=ns_theta, h=h, R=R, rho_in=rho_in, A=A, uniform_disk=uniform_disk, y_tilt=y_tilt, z_tilt=z_tilt, y_tilt2=y_tilt2, z_tilt2=z_tilt2, T_disk = T_disk, spot_disk=spot_disk, T_spot='', spot_beg='', spot_end='', drd_phi=drd_phi, drd_theta=drd_theta, inclination=inclination)
        
elif picture == 1:
    output_filename = 'VIEW.data'



    
#for y_tilt in np.linspace(-10, 40, 90):

#p['y_tilt'] = y_tilt
    
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


#if picture == 0:
#    print 'discostar calculates the lightcurve ...'
#elif picture == 1:
#    print 'discostar draws the picture ...'
    

f = open('./'+direct+'/'+output_filename, "w")    
f_logs = open('./'+logs, "a")

subprocess.call(arg, stdout=f, shell=True)
#subprocess.call(arg, shell=True)

f_logs.write(output_filename+' '+arg+'\n')

f.close
f_logs.close

#print 'done'

if picture == 0:

    print 'lightcurve written to', output_filename
#    print '\n'
#    print '\n'
#    if (spot_disk != 0 and isotrope == 0):
#        print '"../LC2/%s" @legend "%3.0f (%3.0f %2.0f %2.0f) (%3.0f %2.0f %2.0f) h=%1.3f R=%1.2f (%1.1f %1.1f) Td=%2.0f SPOT %d",\\' % (output_filename,z_tilt,PSI_pr,kappa,ns_theta,z_tilt2,y_tilt,y_tilt2,h,R,Lx/1.0e+37,Lx_iso/1.0e+37,T_disk/1000.0,spot_disk)
#    elif (spot_disk != 0 and isotrope == 1):
#        print '"../LC2/%s" @legend "%3.0f ISO (%3.0f %2.0f %2.0f) h=%1.3f R=%1.2f (%1.1f %1.1f) Td=%2.0f SPOT %d",\\' % (output_filename,z_tilt,z_tilt2,y_tilt,y_tilt2,h,R,Lx,Lx_disk/1.0e+37,Lx_iso/1.0e+37,T_disk/1000.0,spot_disk)
#    elif (spot_disk == 0 and isotrope == 1):
#        print '"../LC2/%s" @legend "%3.0f ISO (%3.0f %2.0f %2.0f) h=%1.3f R=%1.2f (%1.1f %1.1f) Td=%2.0f",\\' % (output_filename,z_tilt,z_tilt2,y_tilt,y_tilt2,h,R,Lx/1.0e+37,Lx_iso/1.0e+37,T_disk/1000.0)
#    elif (spot_disk == 0 and isotrope == 0):
#        print '"../LC2/%s" @legend "%3.0f (%3.0f %2.0f %2.0f) (%3.0f %2.0f %2.0f) h=%1.3f R=%1.2f (%1.1f %1.1f) Td=%2.0f",\\' % (output_filename,z_tilt,PSI_pr,kappa,ns_theta,z_tilt2,y_tilt,y_tilt2,h,R,Lx/1.0e+37,Lx_iso/1.0e+37,T_disk/1000.0)   
#    print '\n'
#    print '\n'
#            
        
elif picture == 1:
    print 'data for picture written to', output_filename
