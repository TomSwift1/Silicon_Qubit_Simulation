from qutip import *
import numpy as np
import copy
from numpy import *
import matplotlib.pyplot as plt
import qutip.logging_utils as logging
from pylab import *

# QuTiP control modules
import qutip.control.pulseoptim as cpo


# driven by Hamiltonians- drift (constant) is in z and control (time dependent) is in x


#create state vector corresponding to spin qubit
up = basis(2,0) # i.e. the column vector (1,0)
down = basis(2,1) # i.e. the column vector (0,1)
hbar = 1.054e-34

def H(delta,phi,omega = 42000000):
    """

    Args:
        delta: detuning
        phi: angle of rotation wrt to sigmax and sigmay / radians
        omega: frequency / Hz

    Returns:

    """
    Hamiltonian = -(delta * sigmaz())/2 - (omega/4 * ((np.cos(phi) * sigmax()) - (np.sin(phi) * sigmay())))

    return Hamiltonian

#X GATE
#rotation of pi around x axis
theta_x = np.pi
delta_x = 0
phi_x = 0
w = 42000000 #omega same for y gate
time = theta_x * 2 / w #time same for y gate
time_array = np.linspace(0,time,10000)
print(time_array)

xgate = H(delta_x, phi_x)

#Y GATE
#rotation of pi around y axis
theta_y = np.pi
delta_y = 0
phi_y = np.pi/2

ygate = H(delta_y, phi_y)

#final states under evolution
final_x = mesolve(xgate, up, time_array)
final_y = mesolve(ygate, up, time_array)

#all states of evolution
print(final_x.states)
print(final_y.states)

#final state
print(final_x.states[-1])
print(final_y.states[-1])

#function to do example above automatically
def evolution_H(initial, delta, phi, theta,  omega = 42000000, timestep = 1000):
    """

    Args:
        initial: Qobj. Initial state / ket
        delta: detuning
        phi: rotation axis / radians (pi/2 for y and 0 for x axis)
        omega: angle of rotation
        timestep: desired time step to calculate the time evolution / s

    Returns:

    """

    Hamiltonian = H(delta,phi,omega)
    time = theta * 2 / omega

    time_array = np.linspace(0, time, timestep)
    results = mesolve(Hamiltonian, initial, time_array)
    final_state = results.states[-1]

    return final_state


#FIDELITIES

#define the mixed state
test_state_x = basis(2,1)
initial = basis(2,0)


fidelity_x = fidelity(test_state_x * test_state_x.dag(), evolution_H(initial,0,0,np.pi) * evolution_H(initial,0,0, np.pi).dag())
print("Fidelity of X gate: ", fidelity_x)

test_state_y= -1j*basis(2,1)
fidelity_y = fidelity(test_state_y * test_state_y.dag(), evolution_H(initial,0,np.pi/2, np.pi) * evolution_H(initial,0,np.pi/2, np.pi).dag())
print("Fidelity of Y gate:", fidelity_y)


#DRIFT due to our PERMANENT MAGNET
g = 2
mu = 1.760 * (10)**11 #gyromagnetic ratio electron
B1 = 3e-3 #teslas
delta = g * mu * B1
initial = basis(2,0)

#XGATE
final_state_x = evolution_H(initial,delta,0, np.pi)
print("Fidelity of X gate with a drift:", fidelity(test_state_x * test_state_x.dag(),final_state_x*final_state_x.dag()))


#YGATE FIDELITY ERROR PERMANENT MAGNET
final_state_y = evolution_H(initial,delta,np.pi/2, np.pi)
print("Fidelity of Y gate with a drift:", fidelity(test_state_y * test_state_y.dag(), final_state_y * final_state_y.dag()))


def grid(ind_sep):
    x = np.arange(-2 * ind_sep, 3 * ind_sep, ind_sep)
    xx, yy = np.meshgrid(x, x)

    mat = np.zeros([len(xx), len(xx)])

    for i in range(len(xx)):
        for j in range(len(xx)):
            mat[i, j] = (xx[i][j] ** 2 + yy[i][j] ** 2)

    return mat


sep = 0.2
tot_field_variation = 5*(10**(-3))
final_grid = grid(sep)

b_field = final_grid * tot_field_variation
delta_Y= b_field * 1.7588199959*(10**11)


with np.nditer(delta_Y, op_flags=['readwrite']) as it:
    for x in it:
           x[...] = fidelity(test_state_y * test_state_y.dag(), evolution_H(initial,x,np.pi/2, np.pi) * evolution_H(initial,x,np.pi/2, np.pi).dag())


im = imshow((1-delta_Y)*100,cmap=cm.RdBu,extent=(-2.5*sep, 2.5*sep, 2.5*sep, -2.5*sep))
colorbar(im)
plt.xlabel('x/mm')
plt.ylabel('y/mm')
plt.title('Y-Gate %Error Fidelity')
plt.savefig("Y_gate_fidelity2.png")
plt.show()

#YGATE FIDELITY ERROR SOLENOID

#DRIFT due to solenoid
g = 2
mu = 1.760 * (10)**11 #gyromagnetic ratio electron
B1 = 2e-6 #teslas
delta = (g * mu * B1)
initial = basis(2,0)

def grid(ind_sep):
    x = np.arange(-2 * ind_sep, 3 * ind_sep, ind_sep)
    xx, yy = np.meshgrid(x, x)

    mat = np.zeros([len(xx), len(xx)])

    for i in range(len(xx)):
        for j in range(len(xx)):
            mat[i, j] = (xx[i][j] ** 2 + yy[i][j] ** 2)

    return mat


sep = 0.2
tot_field_variation = 5*(10**(-6))
final_grid = grid(sep)

b_field = final_grid * tot_field_variation
delta_Y= b_field * 1.7588199959*(10**11)


with np.nditer(delta_Y, op_flags=['readwrite']) as it:
    for x in it:
           x[...] = fidelity(test_state_y * test_state_y.dag(), evolution_H(initial,x,np.pi/2, np.pi) * evolution_H(initial,x,np.pi/2, np.pi).dag())


im = imshow((1-delta_Y)*100,cmap=cm.RdBu,extent=(-2.5*sep, 2.5*sep, 2.5*sep, -2.5*sep))
colorbar(im)
plt.xlabel('x/mm')
plt.ylabel('y/mm')
plt.title('Y-Gate %Error Fidelity')
plt.savefig("Y_gate_fidelity_solenoid2.png")
plt.show()