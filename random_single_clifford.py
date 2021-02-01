from qutip.operators import qeye
from qutip.qip.gates import rx,ry
from numpy import pi
import random

def random_single_clifford():
    cliff_gates = [qeye(2),rx(pi/2),ry(pi/2),rx(pi),ry(pi)]
    rand_cliff = random.choice(cliff_gates)
    
    return rand_cliff