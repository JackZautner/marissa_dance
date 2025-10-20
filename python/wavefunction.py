import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Constants (set arbitrary consts to 1 for simplicity)
hbar = 1
m = 1
E = 1
V_0 = 3
k = np.sqrt((2*m*E)/np.power(hbar, 2))
lamb = np.sqrt((2*m*(V_0-E))/np.power(hbar, 2))

# Reflection, transmission coefficients
R = (1-(1j*lamb/k))/(1+(1j*lamb/k))
T = 2/(1+(1j*lamb/k))

def psi_free(x, t):
	return np.exp(1j*k*x-(1j*E*t)/hbar)+R*np.exp(-1j*k*x-(1j*E*t)/hbar)

def psi_trans(x, t):
	return T*np.exp(-lamb*x-(1j*E*t)/hbar)

def psi(x, t):
	psi_full = np.zeros_like(x, dtype=complex)
	psi_full[x<0] = psi_free(x[x<0], t)
	psi_full[x>=0] = psi_trans(x[x>=0], t)
	return psi_full

def step_potential(x):
	V = np.zeros_like(x)
	V = np.split(V, 2)
	V[0][:] = 0
	V[1][:] = V_0
	V = np.array(V).flatten()
	return V

def animate(i):
	plt.clf()
	plt.plot(x, np.real(psi(x, t[i])), label='Real part')
	plt.plot(x, np.imag(psi(x, t[i])), label='Imaginary part')
	plt.plot(x, step_potential(x), label='Step potential')
	plt.xlim(-15, 15)
	plt.ylim(-5, 5)
	plt.legend(loc='upper right')
	plt.xlabel('x')
	plt.ylabel('Wavefunction')

x = np.linspace(-15, 15, 1000)
t = np.linspace(0, 10, 100)

fig = plt.figure()
anm = animation.FuncAnimation(fig, animate, frames=len(t), interval=100)
plt.show()