import numpy as np
import re
#from scipy.interpolate import interp1d
import subprocess
import discostar3_lib as ds3_l
#import matplotlib.pyplot as plt
from scipy.optimize import brute
from scipy.optimize import fmin


data_file = '25_30.dat'

path_to_data_file = './OBSERVED_DATA/'+data_file

border = 0.15

ranges = (slice(3.0e+37, 5.5e+37, 0.5e+37), slice(0.0, 20.0, 2.0), slice(0.0, 30.0, 10.0)) # Lx, y_tilt2, z_tilt2

ranges_disk = (slice(0.0, 20.0, 2.0), slice(0.22, 0.32, 0.02), slice(10000.0, 18000.0, 1000.0)) # y_tilt, R, T_disk

result_disk = brute(ds3_l.residuals_disk, ranges_disk, args=(border, path_to_data_file), full_output=True, finish=None)

result = brute(ds3_l.residuals, ranges, args=(result_disk[0], path_to_data_file,), full_output=True, finish=None)

with open('./temp/optimize_logs', 'a') as logs:
    logs.write(' DATA %s\n Lx\t %e\n y_tilt2\t %f\n z_tilt2\t %f\n y_tilt\t %f\n R\t %f\n T_disk\t %f\n Chisq\t %f\n' % (re.split('\.', data_file)[0], result[0][0], result[0][1], result[0][2], result_disk[0][0], result_disk[0][1], result_disk[0][2], result[1]))
