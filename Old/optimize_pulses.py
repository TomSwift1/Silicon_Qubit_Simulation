import numpy as np
import matplotlib.pyplot as plt
import datetime
from qutip import Qobj, identity, sigmax, sigmaz,sigmay
from qutip.qip import hadamard_transform
import qutip.logging_utils as logging
logger = logging.get_logger()
#Set this to None or logging.WARN for 'quiet' execution
log_level = logging.INFO
#QuTiP control modules
import qutip.control.pulseoptim as cpo

def optimizer(H_d,H_c,U_targ,n_ts,evo_time):
    # H_d = drift Hamiltonian
    # H_c = control Hamiltonian
    # U_targ = target Hamiltonian
    # n_ts = Number of time slots
    # evo_time = Evolution time
    example_name = 'Test'
    
    U_0 = identity(2)
    # Fidelity error target
    fid_err_targ = 1e-10
    # Maximum iterations for the optisation algorithm
    max_iter = 200
    # Maximum (elapsed) time allowed in seconds
    max_wall_time = 120
    # Minimum gradient (sum of gradients squared)
    # as this tends to 0 -> local minima has been found
    min_grad = 1e-20
    
    # pulse type alternatives: RND|ZERO|LIN|SINE|SQUARE|SAW|TRIANGLE|
    p_type = 'RND'
    
    f_ext = "{}_n_ts{}_ptype{}.txt".format(example_name, n_ts, p_type)
    
    result = cpo.optimize_pulse_unitary(H_d, H_c, U_0, U_targ, n_ts, evo_time, 
                fid_err_targ=fid_err_targ, min_grad=min_grad, 
                max_iter=max_iter, max_wall_time=max_wall_time, 
                out_file_ext=f_ext, init_pulse_type=p_type, 
                log_level=log_level, gen_stats=True)
    
    return result.time, result.final_amps