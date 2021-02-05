from qutip import *
import numpy as np
import copy
from numpy import *
import matplotlib.pyplot as plt
import qutip.logging_utils as logging
logger = logging.get_logger()
# Set this to None or logging.WARN for 'quiet' execution
log_level = logging.INFO
# QuTiP control modules
import qutip.control.pulseoptim as cpo


# driven by Hamiltonians- drift (constant) is in z and control (time dependent) is in x


#create state vector corresponding to spin qubit
up = basis(2,0) # i.e. the column vector (1,0)
down = basis(2,1) # i.e. the column vector (0,1)

def H(delta,phi,omega = 42000000):
    """

    Args:
        delta: detuning
        phi: angle of rotation wrt to sigmax and sigmay / radians
        omega: frequency / Hz

    Returns:

    """
    Hamiltonian = -(delta * sigmaz())/2 - omega/2 * (np.cos(phi) * sigmax() - np.sin(phi) * sigmay())

    return Hamiltonian

#X GATE
#rotation of pi around x axis
theta_x = np.pi
delta_x = 0
phi_x = 0
w = 65000000 #omega same for y gate
time = theta_x * 2 / w #time same for y gate
time_array = np.linspace(0,time,10000)
print(time_array)

xgate = H(delta_x, phi_x)

#Y GATE
#rotation of pi around y axis
theta_y = np.pi
delta_y = 0
phi_y = 0

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
def evolution_H(initial, delta, phi,  omega = 42000000, timestep = 1000):
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
    time = theta_x * 2 / omega

    time_array = np.linspace(0, time, timestep)
    results = mesolve(Hamiltonian, initial, time_array)
    final_state = results.states[-1]

    return final_state


#FIDELITIES

#define the mixed state
plus_state = (basis(2,0)+basis(2,1))/np.sqrt(2)
fidelity_x = fidelity(plus_state * plus_state.dag(), evolution_H(up,0,0) * evolution_H(up,0,0).dag())
print(fidelity_x)

fidelity_y = fidelity(plus_state * plus_state.dag(), evolution_H(up,0,np.pi/2) * evolution_H(up,0,np.pi/2).dag())
print(fidelity_y)


#DRIFT due to our PERMANENT MAGNET
g = 2
mu = 1.760 * (10)**11 #gyromagnetic ratio electron
delta = g * mu * deltaB
initial = basis(2,0)

Hamiltonian = H(delta, np.pi)

