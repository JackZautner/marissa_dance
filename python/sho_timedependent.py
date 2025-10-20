import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# CONSTANTS
hbar = 1

# State info
m = 1
a = 1
d = 1
omega = 1

# Physics
def potential(x, t):
	x = np.complex128(x)
	return 0.5*m*omega**2*x**2

def prefactor(x, t):
	x = np.complex128(x)
	gauss_factor = 1 / (np.power(np.pi, 0.25) * np.sqrt(d))
	oscillating_factor = 1 / np.sqrt(np.cos(omega*t)+(1j*hbar)/(m*omega*d**2)*np.sin(omega*t))
	return gauss_factor*oscillating_factor

def exponential_one(x, t):
	x = np.complex128(x)
	argument_1 = (1j*m*omega*x**2)/(2*hbar*np.tan(omega*t))
	argument_2 = (1j*m*omega*x*a)/(hbar*np.sin(omega*t))
	argument_3 = (1j*m*omega*a**2)/(2*hbar*np.tan(omega*t))
	argument = ( argument_1 - argument_2 + argument_3 )
	return np.exp(argument)

def exponential_two(x, t):
	x = np.complex128(x)
	argument = - ( (m*omega*x)/(hbar*np.sin(omega*t)) - (a*m*omega)/(hbar*np.tan(omega*t)) )**2 / 2 * ( 1/d**2 - (1j*m*omega)/(hbar*np.tan(omega*t)) )
	return np.exp(argument)

def psi(x, t):
	x = np.complex128(x)
	return prefactor(x, t) * exponential_one(x, t) * exponential_two(x, t)

def probability(x, t):
	return np.conjugate(psi(x,t))*psi(x,t)

def test_wave(x, t):
	return np.exp( 1j*x - 1j*omega*t)

# Animation stuff
def animate(i):
	print(i)
	plt.clf()

	plt.plot(x, np.real(psi(x, t[i])), label='Real part')
	plt.plot(x, np.imag(psi(x, t[i])), label='Imaginary part')
	#plt.plot(x, probability(x, t[i]), label='Probability distribution')
	plt.plot(x, potential(x, t[i]), label='Oscillator potential')

	plt.xlim(-5, 5)
	plt.ylim(-2, 2)

	plt.legend(loc='upper right')
	plt.xlabel('x')
	plt.ylabel('psi(x, t)')

# Coordinates
x = np.linspace(-5, 5, 10000)
t = np.linspace(1, 10, 500)

fig = plt.figure()
anm = animation.FuncAnimation(fig, animate, frames=len(t), interval=100)
plt.show()