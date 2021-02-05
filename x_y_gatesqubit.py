from qutip import *
from qutip.qip.device import Processor
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

#X gate - target H is sigma x

# processor == quantum device
# driven by Hamiltonians- drift (constant) is in z and control (time dependent) is in x



#create state vector corresponding to spin qubit
up = basis(2,0) # i.e. the column vector (1,0)
down = basis(2,1) # i.e. the column vector (0,1)

#define a function of rotations
def rot(theta, direction):
    """

    Args:
        omega: frequency calculated (ask Tom)/ Hz
        time: time driving the Hamiltonian needed to obtain the desired gate / s
        direction: direction of rotation in the Bloch sphere. 'x' or 'y'

    Returns:
        QuTip Obj: matrix of rotation

    """

    if direction == 'x':
        rotation = Qobj(
            [
                [
                    np.cos(theta/2),
                    -np.sin(theta/2 )*1j,
                ],
                [
                    -np.sin(theta/2 )*1j,
                    np.cos(theta/2 )],
            ]
        )
        return rotation

    if direction == 'y':
        rotation = Qobj(
            [
                [
                    np.cos((omega * time)/2 ),
                    - np.sin((omega * time)/2 )],
                [
                    np.sin((omega * time)/2 ),
                    np.cos((omega * time)/2 )],

            ]
        )
        return rotation



rotation_x = rot(np.pi/2,'x') #QuTip returns not Hermitian ????
print(rotation_x * rotation_x.dag())
print(rotation_x)

#X GATE
#theta is pi
#look at time calculations
xgate = rot(1e-07, 'x')

#YGATE
#theta is pi
#look at time calculations
ygate = rot(1e-07, 'y')

#xgate and ygate have both a time of 10^-7
#refer to calculations for this

print(xgate * up)
print(xgate * down)

print(ygate * up)
print(ygate * down)

def density_matrix (operator, state):
    """

    Args:
        operator: QObj. Matrix of operator
        state: QObj. Ket of initial state

    Returns: QObj. Density matrix of final state

    """
    finalstate = operator * state

    dens = finalstate * finalstate.dag()

    return dens

#trial for fidelities

dens1 = density_matrix(xgate,up)
dens2 = density_matrix(xgate,down)

print(fidelity(dens1,dens1))