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

# processor == quantum device
# driven by Hamiltonians- drift (constant) is in z and control (time dependent) is in x

# Constant Hamiltonian
H_d = sigmaz()
# The control Hamiltonian
H_c = [sigmax()]
# start point for the gate evolution
U_0 = qeye(2)
# random sequence of clifford gates
# repeated for the following lengths
# 1,12,23,34,45,56,67,78,89,100
U_targ = random_clifford_sequence(1)