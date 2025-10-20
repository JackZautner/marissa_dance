import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import warnings
warnings.filterwarnings('ignore')

'''
PHYSICAL CONSTANTS
'''
nu_0 = 2.584 # Resonant wavenumber
dnu = 0.0333 # Discretization of wavenumber
V_0 = 10 # Potential
hbar = 1 # Reduced Planck constant
m = 0.5 # Mass
L = 5 # Length of well
a = np.power(8, -0.5) # Standard deviation of wavepacket peak

j = np.arange(-100, 101) # Index range for wavenumbers
nu = nu_0 + j * dnu # Wavenumbers

'''
GENERIC METHODS, CONSTANTS
'''
xticks = 1000
xlim_up = 35

def val_to_index(value, ticks, lim, floor=True):
	if floor:
		return value*np.floor_divide(ticks, lim)
	else:
		return value*int(np.ceil(ticks / lim))

def square_mod(fn):
	return np.real(fn*np.conjugate(fn))

'''
SCATTERING PHYSICS
'''
def k(nu):
	return np.sqrt(np.power(nu, 2)+V_0)

def reflection_coefficient(nu):
	return np.exp(-2j*L*nu)*(1j*nu*np.tan(k(nu)*L)+k(nu))/(1j*nu*np.tan(k(nu)*L)-k(nu))

def well_amplitude(nu):
	return (-2j*nu*np.exp(-1j*L*nu))/(k(nu)*np.cos(k(nu)*L)-1j*nu*np.sin(k(nu)*L))

def gauss(nu):
	# Gaussian envelope for wavepacket
	return np.exp(-np.power(nu-nu_0, 2)/(2*np.power(a, 2)))

def potential(x):
	V = np.zeros_like(x)
	we = np.floor_divide(1000, 7)
	V[0:we] = -1 # Well
	V[we:] = 0 # Potential region
	return V

def psi_out(x, t):
	psi = None
	for nu_j in nu:
		R = reflection_coefficient(nu_j)
		psi_j = dnu*gauss(nu_j)*(np.exp(-1j*nu_j*x)+R*np.exp(1j*nu_j*x))*np.exp(-(1j*hbar*np.power(nu_j, 2)*t)/(2*m))
		if psi is None:
			psi = psi_j
		else:
			psi += psi_j
	return psi

def psi_in(x, t):
	psi = None
	for nu_j in nu:
		A = well_amplitude(nu_j)
		psi_j = dnu*gauss(nu_j)*A*np.sin(np.sqrt(np.power(nu_j, 2)+V_0)*x)*np.exp(-(1j*hbar*np.power(nu_j, 2)*t)/(2*m))
		if psi is None:
			psi = psi_j
		else:
			psi += psi_j
	return psi

'''
ANIMATION
'''
x = np.linspace(0, xlim_up, xticks) # x-axis
t = np.linspace(-3, 3, 100) # t-axis

adj_L = val_to_index(L, xticks, xlim_up) # Index-adjusted L, loored version

x_in = x[0:adj_L+1]
x_out = x[adj_L:]

print(psi_in(x, 0)[adj_L], np.conjugate(psi_in(x, 0)[adj_L]))

def animate(i):
	plt.clf()
	pdf_in = square_mod(psi_in(x, t[i]))
	pdf_out = square_mod(psi_out(x, t[i]))
	plt.plot(x_out, pdf_out[adj_L:], label='Probability', color='blue')
	plt.plot(x_in, pdf_in[0:adj_L+1], color='blue')
	print(f'in: {np.gradient(pdf_in)[adj_L+1]}, out: {np.gradient(pdf_out)[adj_L]}')
	plt.plot(x, potential(x), label='Potential', color='red')
	plt.xlim(0, 35)
	plt.ylim(-1, 2)
	plt.legend(loc='upper right')
	plt.xlabel('x')
	plt.ylabel('Probability')

fig = plt.figure()
anm = animation.FuncAnimation(fig, animate, frames=len(t), interval=100, repeat=False)
plt.show()

anm.save('.\\animation_2.0.gif', writer='pillow')