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

def H(delta,phi,omega = 65000000):
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
time_array = np.linspace(0,time,1000)
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

