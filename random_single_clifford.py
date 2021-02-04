from qutip.operators import qeye
from qutip.qip.operations import rx,ry
from numpy import pi
import random
from qutip.states import basis
from numpy import pi


def random_single_clifford():
    cliff_gates = [qeye(2),rx(pi/2),ry(pi/2),rx(pi),ry(pi)]
    rand_cliff = random.choice(cliff_gates)
    
    return rand_cliff


def random_clifford_sequence(n_gates):
    
    sequence = [None]*n_gates
    final = qeye(2)
    for i in range(n_gates-1):
        a_random_clifford = random_single_clifford()
        sequence[i] = a_random_clifford
        final = a_random_clifford*final

    sequence[n_gates-1] = final.inv()
    
    return sequence