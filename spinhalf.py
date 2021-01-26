from qutip import *
import matplotlib.pyplot as plt
from numpy import sqrt, pi, array, sin, cos, arange

#convention p means positive and m means negative
pz = basis(2,0) # i.e. the column vector (1,0)
mz = basis(2,1) # i.e. the column vector (0,1)
px = 1/sqrt(2)*(pz + mz)
mx = 1/sqrt(2)*(pz - mz)
py = 1/sqrt(2)*(pz + 1j*mz)
my = 1/sqrt(2)*(pz - 1j*mz)

print(pz.dag() * mz) #as expected this is 0

#spin 1/2 projection operators
Sx = 1/2.0*sigmax()
Sy = 1/2.0*sigmay()
Sz = 1/2.0*sigmaz()

#Hamiltonian in a B field
gamma = 1
B = 2
H = -gamma * Sz * B #if B is a classical field there is another way to do this

#if B is a classical field
omega = 5
Hz = -omega*Sz

#evolution of the system
t = arange(0,4*pi/omega,0.05)
expect_ops = [Sx,Sy,Sz]
psi0 = px
result = sesolve(Hz, psi0, t, expect_ops)

labels = ['x', 'y', 'z']
style = {'x':'-.', 'y':'--', 'z':'-'}
for r,l in zip(result.expect,labels):
    plt.plot(t*omega/pi, r, style[l], label="$\langle S_%c \\rangle $" % l)
    plt.xlabel("Time ($\Omega t/\pi$)", size=18)
    plt.legend(fontsize=16)
plt.show()

