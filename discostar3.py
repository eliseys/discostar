import numpy as np
#from scipy.interpolate import interp1d
import subprocess
import discostar3_lib as ds3_l
import matplotlib.pyplot as plt
from scipy.optimize import brute
from scipy.optimize import fmin


data_file = "./OBSERVED_DATA/25_30.dat"

with open(data_file, 'r') as data:
    x_obs, y_obs = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)

def residuals(pars):

    y_tilt = pars[0]
    y_tilt2 = pars[1]
    z_tilt = pars[2]
    
    f = ds3_l.light_curve(y_tilt, y_tilt2, z_tilt)
    r = (y_obs - f(x_obs))**2

    return r.sum()


ranges = (slice(0.0, 30.0, 5.0), slice(0.0, 30.0, 5.0), slice(0.0, 50.0, 10.0))
result = brute(residuals, ranges, full_output=True)

print "RESULT", result[0], result[1]

args = result[0]

#f = ds3_l.light_curve(y_tilt, y_tilt2, z_tilt)
f = ds3_l.light_curve(*args)

x = np.arange(0, 1, 0.01)
y = f(x)

plt.plot(x, y, '-', x_obs, y_obs, '.')
#plt.show()
plt.savefig('./temp/phases_27_10_17/20_25.png')




