import numpy as np
from scipy.fft import fft
import matplotlib.pyplot as plt

def gaussian(x, A, mu, sigma):
	'''Gaussian distribution (unnormalized)'''
	return A*np.exp( -(x - mu)**2 / (2*sigma**2) )

x_axis = np.linspace(-10, 10, 1000)
fn = gaussian(x_axis, 1, 0, 1)
chi_fn = fft(fn)

plt.plot(x_axis, chi_fn)
plt.show()