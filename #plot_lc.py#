import numpy as np
import matplotlib.pyplot as plt

with open("TEST.data", "r") as f:
    phase, flux = np.loadtxt(f, dtype=('float'), usecols=(0,1), unpack=True)



s = len(flux)


#phase = np.roll(phase, s//2)
flux = np.roll(flux, s//2)



print(phase)
print(flux)
    
plt.plot(phase, flux)
plt.grid()
plt.show()
