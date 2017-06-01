import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import basinhopping
#from scipy.stats import chisquare
from scipy.optimize import brute
#from scipy.optimize import fmin_l_bfgs_b
#from scipy.optimize import fmin_powell
from scipy.optimize import fmin
import matplotlib.pyplot as plt
import math
import subprocess
import opt_lib

# initial guess values in the parameter file
inguess_LT = opt_lib.inguess_LT
inguess_yz = opt_lib.inguess_yz
result_LT = np.array([1.0, 1.0])
result_yz = np.array([1.0, 1.0, 1.0])

width = 0.2
borders_LT = [((1 - width),(1 + width)),((1 - width),(1 + width))]
borders_yz = [((1 - width),(1 + width)),((1 - width),(1 + width)),((1 - width),(1 + width))]

for z_tilt in range(0,18,18):

    with open("output.data", 'a') as out:
        out.write("PHASE\t{data_files}\t{z_tilt}\n".format(data_files = opt_lib.data_files[z_tilt], z_tilt = z_tilt))    

    with open(opt_lib.data_directory+opt_lib.data_files[z_tilt], 'r') as data:
        opt_lib.x_obs, opt_lib.y_obs = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)
        
        for i in range(30):
            #result_LT = brute(opt_lib.fit_LT, borders_LT, args=(result_yz[0] * inguess_yz[0], result_yz[1] * inguess_yz[1], result_yz[2] * inguess_yz[2], z_tilt), Ns=3)

            #r_2_sum = opt_lib.fit_LT(result_LT, result_yz[0] * inguess_yz[0], result_yz[1] * inguess_yz[1], result_yz[2] * inguess_yz[2], z_tilt)
            
            #with open("output.data", 'a') as out:
            #    out.write("{i}\t LT brute\t{Lx}\t{T_disk}\t{y_tilt}\t{y_tilt2}\t{z_tilt2}\t{r_2_sum}\n".format(i = i,
            #                                                                                                  Lx = result_LT[0] * inguess_LT[0],
            #                                                                                                  T_disk = result_LT[1] * inguess_LT[1],
            #                                                                                                  y_tilt = result_yz[0] * inguess_yz[0],
            #                                                                                                  y_tilt2 = result_yz[1] * inguess_yz[1],
            #                                                                                                  z_tilt2 = result_yz[2] * inguess_yz[2],
            #                                                                                                  r_2_sum = r_2_sum))

                
            result_LT = fmin(opt_lib.fit_LT, result_LT, args=(result_yz[0] * inguess_yz[0], result_yz[1] * inguess_yz[1], result_yz[2] * inguess_yz[2], z_tilt))

            r_2_sum = opt_lib.fit_LT(result_LT, result_yz[0] * inguess_yz[0], result_yz[1] * inguess_yz[1], result_yz[2] * inguess_yz[2], z_tilt)

            with open("output.data", 'a') as out:
                out.write("{i}\t LT fmin\t{Lx}\t{T_disk}\t{y_tilt}\t{y_tilt2}\t{z_tilt2}\t{r_2_sum}\n".format(i = i,
                                                                                                              Lx = result_LT[0] * inguess_LT[0],
                                                                                                              T_disk = result_LT[1] * inguess_LT[1],
                                                                                                              y_tilt = result_yz[0] * inguess_yz[0],
                                                                                                              y_tilt2 = result_yz[1] * inguess_yz[1],
                                                                                                              z_tilt2 = result_yz[2] * inguess_yz[2],
                                                                                                              r_2_sum = r_2_sum))


                
            #result_yz = brute(opt_lib.fit_yz, borders_yz, args=(result_LT[0] * inguess_LT[0], result_LT[1] * inguess_LT[1], z_tilt), Ns=3)
     
            #r_2_sum = opt_lib.fit_yz(result_yz, result_LT[0] * inguess_LT[0], result_LT[1] * inguess_LT[1], z_tilt)

            #with open("output.data", 'a') as out:
            #    out.write("{i}\t yz brute\t{Lx}\t{T_disk}\t{y_tilt}\t{y_tilt2}\t{z_tilt2}\t{r_2_sum}\n".format(i = i,
            #                                                                                                  Lx = result_LT[0] * inguess_LT[0],
            #                                                                                                  T_disk = result_LT[1] * inguess_LT[1],
            #                                                                                                  y_tilt = result_yz[0] * inguess_yz[0],
            #                                                                                                  y_tilt2 = result_yz[1] * inguess_yz[1],
            #                                                                                                  z_tilt2 = result_yz[2] * inguess_yz[2],
            #                                                                                                  r_2_sum = r_2_sum))


            result_yz = fmin(opt_lib.fit_yz, result_yz, args=(result_LT[0] * inguess_LT[0], result_LT[1] * inguess_LT[1], z_tilt))

            r_2_sum = opt_lib.fit_yz(result_yz, result_LT[0] * inguess_LT[0], result_LT[1] * inguess_LT[1], z_tilt)

            with open("output.data", 'a') as out:
                out.write("{i}\t yz fmin\t{Lx}\t{T_disk}\t{y_tilt}\t{y_tilt2}\t{z_tilt2}\t{r_2_sum}\n".format(i = i,
                                                                                                             Lx = result_LT[0] * inguess_LT[0],
                                                                                                             T_disk = result_LT[1] * inguess_LT[1],
                                                                                                             y_tilt = result_yz[0] * inguess_yz[0],
                                                                                                             y_tilt2 = result_yz[1] * inguess_yz[1],
                                                                                                             z_tilt2 = result_yz[2] * inguess_yz[2],
                                                                                                             r_2_sum = r_2_sum))


                
            #******************************************************************
            
            #f = flx(Lx, y_tilt, T_disk, y_tilt2, z_tilt2, z_tilt)
            #xnew = np.linspace(0.0, 1.0, num=500, endpoint=True)
            #plt.plot(x_obs, y_obs, '.', xnew, f(xnew), '-')
            #plt.show()
        
