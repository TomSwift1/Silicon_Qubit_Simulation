from numpy import cos,sin,pi,sqrt
from qutip import *
import numpy as np

hbar = 1
w_1 = 60*10**6

def spin_hamiltonian(delta,phi):
    
    a = -(hbar/2)*delta*sigmaz()
    b = -(hbar*w_1/4)*cos(phi)*sigmax()
    c = (hbar*w_1/4)*sin(phi)*sigmay()
    
    H = a + b + c
    
    return H

def hadamard_evolution(delta):
    plus_state = (basis(2,0)+basis(2,1))/sqrt(2)
    
    #Ry by pi/2 rotation
    psi = basis(2,0)
    phi = pi/2
    theta = pi/2
    state1 = general_evolution(psi,theta,phi,delta)
    
    #Rx by pi rotation
    phi = 0
    theta = pi
    final_state = general_evolution(state1,theta,phi,delta)
    
    fid = fidelity(plus_state*plus_state.dag(),final_state*final_state.dag())
    
    return fid

def x_evolution(delta):
    test_state = -1j*basis(2,1)
    
    #Rx by pi rotation
    psi = basis(2,0)
    phi = 0
    theta = pi
    final_state = general_evolution(psi,theta,phi,delta)
    
    fid = fidelity(test_state*test_state.dag(),final_state*final_state.dag())
    
    return fid

def opt_x_evolution(delta):
    test_state = (1j*basis(2,1) + basis(2,0))/sqrt(2)
    #Following Jonathon Jones Optimised 90x
    
    #Rx by 385 rotation
    psi = basis(2,0)
    phi = 0
    theta = pi*(385/180)
    state1 = general_evolution(psi,theta,phi,delta)

    #R -x by pi rotation
    phi = pi
    theta = pi*(320/180)
    state2 = general_evolution(state1,theta,phi,delta)
    
    #Rx by pi rotation
    phi = 0
    theta = pi*(25/180)
    final_state = general_evolution(state2,theta,phi,delta)
    
    fid = fidelity(test_state*test_state.dag(),final_state*final_state.dag())
    
    return fid

def z_evolution(delta):
    
    test_state = 1j*basis(2,0)
    
    #R(-X) by pi/2 rotation
    psi = basis(2,0)
    phi = pi
    theta = pi/2
    state1 = general_evolution(psi,theta,phi,delta)
    
    #RY by pi rotation
    phi = pi/2
    theta = pi
    state2 = general_evolution(state1,theta,phi,delta)
    
    #Rx by pi/2 rotation
    phi = 0
    theta = pi/2
    final_state = general_evolution(state2,theta,phi,delta)
    
    fid = fidelity(test_state*test_state.dag(),final_state*final_state.dag())
    
    return fid

def general_evolution(psi,theta,phi,delta):
    
    ev_time = 2*theta/w_1
    H = spin_hamiltonian(delta,phi)
    t = np.linspace(0,ev_time,10)
    
    result = mesolve(H,psi,t)
    final_state = result.states[-1]
    
    return final_state

def create_grid(sep):
    x = np.arange(-2*sep,3*sep,sep)
    xx,yy = np.meshgrid(x,x)
    
    mat = np.zeros([len(xx),len(xx)])
    
    for i in range(len(xx)):
        for j in range(len(xx)):
            mat[i,j] = (xx[i][j]**2 + yy[i][j]**2)
            
    return mat
    