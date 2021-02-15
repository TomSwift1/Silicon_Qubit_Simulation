from qutip import *
import numpy as np
import copy
from numpy import *
import matplotlib.pyplot as plt
import qutip.logging_utils as logging
from pylab import *
import seaborn as sns
from X_Y_gates_evolution import H, evolution_H

#creating array of field variations

field_var = np.linspace(1e-10,0.001,100)
length = np.arange(0, len(field_var),1)

g = 2
mu = 1.760 * (10)**11 #gyromagnetic ratio electron
delta_array= []

for i in length:
    delta = g * mu * field_var[i]
    delta_array.append(delta)

print(delta_array)

test_state_x = basis(2,1)
initial = basis(2,0)
fidelity_array = []

for i in length:
    fide = fidelity(test_state_x * test_state_x.dag(),
             evolution_H(initial, delta_array[i], 0, np.pi) * evolution_H(initial, delta_array[i], 0, np.pi).dag())
    fidelity_array.append(fide)

print(fidelity_array)

sns.set_style("whitegrid")

plt.figure()
plt.plot(field_var,fidelity_array)
plt.title("Fidelity Dependence on Field Variation")
plt.xlabel("Field Variation / T/mm")
plt.ylabel("Fidelity")

plt.savefig("FidelityvsField.png")
plt.show()


