import numpy as np
import re
#from scipy.interpolate import interp1d
import subprocess
import discostar3_lib as ds3_l
#import matplotlib.pyplot as plt
from scipy.optimize import brute
from scipy.optimize import fmin


#phase = {
#    '25_30.dat':[0,	270],
#    '20_25.dat':[342,	288],
#    '15_20.dat':[324,	306],
#    '10_15.dat':[306,	324],
#    '05_10.dat':[288,	342],
#    '00_05.dat':[270,	0],
#    '95_00.dat':[252,	18],
#    '90_95.dat':[234,	36],
#    '85_90.dat':[216,	54],
#    '80_85.dat':[198,	72],
#    '75_80.dat':[180,	90],
#    '70_75.dat':[162,	108],
#    '65_70.dat':[144,	126],
#    '60_65.dat':[126,	144],
#    '55_60.dat':[108,	162],
#    '50_55.dat':[90,	180],
#    '45_50.dat':[72,    198],
#    '40_45.dat':[54,	216],
#    '35_40.dat':[36,	234],
#    '30_35.dat':[18,	252]
#}

phase = {
    '30_35.dat':[18,	252]
}


border = 0.15

y_tilt = [10.0, 11.0, 1.0] # min, max, step
#R =      [0.22, 0.24, 0.02]
T_disk = [12000.0, 13000.0, 1000.0]

Lx =      [3.5e+37, 4.0e+37, 0.5e+37]
y_tilt2 = [10.0, 11.0, 1.0]
z_tilt2 = [30.0, 40.0, 10.0]

#ranges_disk = (slice(y_tilt[0], y_tilt[1], y_tilt[2]), slice(R[0], R[1], R[2]), slice(T_disk[0], T_disk[1], T_disk[2])) # y_tilt, R, T_disk
ranges_disk = (slice(y_tilt[0], y_tilt[1], y_tilt[2]), slice(T_disk[0], T_disk[1], T_disk[2])) # y_tilt, T_disk


ranges = (slice(Lx[0], Lx[1], Lx[2]), slice(y_tilt2[0], y_tilt2[1], y_tilt2[2]), slice(z_tilt2[0], z_tilt2[1], z_tilt2[2])) # Lx, y_tilt2, z_tilt2





for i, x in enumerate(phase):

    zpsi = phase[x] # z_tilt & PSI_pr
    data_file = x

    path_to_data_file = './OBSERVED_DATA/'+data_file

    result_disk = brute(ds3_l.residuals_disk, ranges_disk, args=(border, zpsi, path_to_data_file), full_output=True, finish=None)

    result = brute(ds3_l.residuals, ranges, args=(result_disk[0], zpsi, path_to_data_file,), full_output=True, finish=None)

    

    
    with open('./temp/optimize_logs', 'a') as logs:
        logs.write('\nDATA %s\nLx\t\t %e\nPSI_pr\t\t %f\nz_tilt\t\t %f\ny_tilt\t\t %f\ny_tilt2\t\t %f\nz_tilt2\t\t %f\nT_disk\t\t %f\nChisq\t\t %f\n'
                   % (re.split('\.',data_file)[0],
                      result[0][0],
                      zpsi[1],
                      zpsi[0],
                      result_disk[0][0],
                      result[0][1],
                      result[0][2],
                      result_disk[0][1],
                      result[1]))


    print '-------', re.split('\.',data_file)[0], '-------' 




        

