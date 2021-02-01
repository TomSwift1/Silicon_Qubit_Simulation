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

# Constant Hamiltonian
H_d = sigmaz()
# The control Hamiltonian
H_c = [sigmax()]
# start point for the gate evolution
U_0 = qeye(2)
# Target for the gate evolution X gate
U_targ = sigmax()

# Defining time evolution parameters
# Number of time slots
n_ts = 1000
# Time allowed for the evolution
evo_time = 10

# Conditions that will make pulse terminate
# Fidelity error target
fid_err_targ = 1e-10
# Maximum iterations for the optisation algorithm
max_iter = 500
# Maximum (elapsed) time allowed in seconds
max_wall_time = 30
# Minimum gradient (sum of gradients squared)
# as this tends to 0 -> local minima has been found
min_grad = 1e-20

# Initial Pulse Type
# pulse type alternatives: RND|ZERO|LIN|SINE|SQUARE|SAW|TRIANGLE|
p_type = 'SINE'

# Output files
# uncomment to create output files
# f_ext = "{}_n_ts{}_ptype{}.txt".format(example_name, n_ts, p_type)
f_ext = None

# Optimisation
# the results of this is time and initial_amps
result = cpo.optimize_pulse_unitary(H_d, H_c, U_0, U_targ, n_ts, evo_time,
                                    fid_err_targ=fid_err_targ, min_grad=min_grad,
                                    max_iter=max_iter, max_wall_time=max_wall_time,
                                    out_file_ext=f_ext, init_pulse_type=p_type,
                                    log_level=log_level, gen_stats=True)
print(result.time.shape)
print(result.initial_amps.shape)


# plot the results
fig1 = plt.figure()
ax1 = fig1.add_subplot(2, 1, 1)
ax1.set_title("Initial control amps")
ax1.set_xlabel("Time")
ax1.set_ylabel("Control amplitude")
ax1.step(result.time, np.hstack((result.initial_amps[:, 0], result.initial_amps[-1, 0])), where='post')

ax2 = fig1.add_subplot(2, 1, 2)
ax2.set_title("Optimised Control Sequences")
ax2.set_xlabel("Time")
ax2.set_ylabel("Control amplitude")
ax2.step(result.time,
         np.hstack((result.final_amps[:, 0], result.final_amps[-1, 0])),
         where='post')
plt.tight_layout()
plt.savefig("pulseSINE.png")
plt.show()

# PROCESSOR

T1 = 2 #relaxation time
T2 = 10 #dephasing

processor = Processor(N=1, spline_kind="step_func",t1=T1,t2=T2)
processor.add_control(sigmax(), 0)
processor.add_drift(sigmaz(), 0)

plt.figure()
processor.pulses[0].tlist = result.time
processor.pulses[0].coeff = result.initial_amps[:,0]
processor.plot_pulses()
plt.show()
