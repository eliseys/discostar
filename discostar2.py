#from __future__ import print_function
import numpy as np
import subprocess

import gspread
from oauth2client.service_account import ServiceAccountCredentials

direct = "./LC2" # directory for LC`s

scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name('../discostar_parameters/client_secret.json', scope)
client = gspread.authorize(creds)
#sheet = client.open("Discostar").sheet1
sheet = client.open("Discostar").get_worksheet(1)

parameter_name = sheet.row_values(1) # 1st row in google spreadsheet -- names of the parameters

parameter_name = [x for x in parameter_name if x != ''] # deleting empty strings

p = {}

for i, x in enumerate(parameter_name):
    p[x] = sheet.col_values(i+1)[1:]
    p[x] = [u for u in p[x] if u != '']
    p[x] = [u if x == 'DATA' else float(u) for u in p[x]]


#for x in parameter_name:
#    print p[x][0]


for i, x in enumerate(p['DATA']):

    z_tilt = p['z_tilt'][i]
    Lx = p['Lx'][i]
    y_tilt = p['y_tilt'][i]
    y_tilt2 = p['y_tilt2'][i]
    z_tilt2 = p['z_tilt2'][i]
    T_disk = p['T_disk'][i]
    h = p['h'][i]
    R = p['R'][i]
    picture = p['picture'][i]
    inclination = p['inclination'][i]
    q = p['q'][i]
    mu = p['mu'][i]
    beta = p['beta'][i]
    u = p['u'][i]
    albedo = p['albedo'][i]
    lc_num = p['lc_num'][i]
    star_tiles = p['star_tiles'][i]
    disk_tiles = p['disk_tiles'][i]
    threads = p['threads'][i]
    T_star = p['T_star'][i]
    lambda_A = p['lambda_A'][i]
    a = p['a'][i]
    PSI_pr = p['PSI_pr'][i]
    kappa = p['kappa'][i]
    isotrope = p['isotrope'][i]
    Lx_disk = p['Lx_disk'][i]
    spot_disk = p['spot_disk'][i]
    T_spot = p['T_spot'][i]
    spot_beg = p['spot_beg'][i]
    spot_end = p['spot_end'][i]
    ns_theta = p['ns_theta'][i]
    spot_rho_in = p['spot_rho_in'][i]
    spot_rho_out = p['spot_rho_out'][i]
    drd_phi = p['drd_phi'][i]
    drd_theta = p['drd_theta'][i]
    


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
           str(p['q'][i]) + ' ' +
           str(p['mu'][i]) + ' ' +
           str(p['beta'][i]) + ' ' +  
           str(p['u'][i]) + ' ' + 
           str(p['albedo'][i]) + ' ' +
           str(p['Lx'][i]) + ' ' + 
           str(p['h'][i]) + ' ' + 
           str(p['R'][i]) + ' ' + 
           str(p['y_tilt'][i]) + ' ' + 
           str(p['z_tilt'][i]) + ' ' + 
           str(p['picture'][i]) + ' ' +
           str(p['inclination'][i]) + ' ' + 
           str(int(p['lc_num'][i])) + ' ' + 
           str(int(p['star_tiles'][i])) + ' ' + 
           str(int(p['disk_tiles'][i])) + ' ' +
           str(int(p['threads'][i])) + ' ' +
           str(p['T_disk'][i]) + ' ' + 
           str(p['T_star'][i]) + ' ' + 
           str(p['lambda_A'][i]) + ' ' +
           str(p['a'][i]) + ' ' +
           str(p['y_tilt2'][i]) + ' ' + 
           str(p['z_tilt2'][i]) + ' ' +
           str(p['PSI_pr'][i]) + ' ' +
           str(p['kappa'][i]) + ' ' +
           str(p['isotrope'][i]) + ' ' +
           str(p['Lx_disk'][i]) + ' ' +
           str(p['spot_disk'][i]) + ' ' +
           str(p['T_spot'][i]) + ' ' +
           str(p['spot_beg'][i]) + ' ' +
           str(p['spot_end'][i]) + ' ' +
           str(p['ns_theta'][i]) + ' ' +
           str(p['spot_rho_in'][i]) + ' ' +
           str(p['spot_rho_out'][i]) + ' ' +
           str(p['drd_phi'][i]) + ' ' +
           str(p['drd_theta'][i])
    )


    if picture == 0:
        print 'discostar calculates the lightcurve ...'
    elif picture == 1:
        print 'discostar draws the picture ...'
        


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

