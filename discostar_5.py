import lmfit
import numpy as np
import discostar4_lib as d4l

z_tilt0 = 90.0


lc_num = 50


z_tilt_step = 360.0/(20.0*lc_num)


PSI_pr = 0.0
Lx = 1.12e+37
y_tilt0 = 39.0

y_tilt_step = (39.0 - 25.0)/lc_num



y_tilt2 = 38.0
z_tilt2 = -20.0
disk_flux = 0
T_disk = 32000.0
theta = -3.0
kappa = 195.0


for i in range(lc_num + 1): 

    z_tilt = z_tilt0 - z_tilt_step * i
    y_tilt = y_tilt0 - y_tilt_step * i

    #PSI_pr = 0.0 + z_tilt_step*i

                
    f = d4l.lc_disk(z_tilt, PSI_pr, Lx, y_tilt, y_tilt2, z_tilt2, disk_flux, T_disk, theta, kappa)

    #print '{:8f}'.format(np.sum(np.square(data_short - f(x_short)))/N), '\t', '{:8f}'.format(z_tilt), '\t', '{:8f}'.format(y_tilt),  '\t','{:8f}'.format(y_tilt2), '\t', '{:8f}'.format(z_tilt2), '\t', '{:10f}'.format(T_disk), '\t', '{:E}'.format(Lx)

    Flux = float(f(float(i)/lc_num))
    
    
    print '{:8f}'.format(float(i)/lc_num - 0.5), '\t', '{:8f}'.format(Flux)

