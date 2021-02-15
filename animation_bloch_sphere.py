import matplotlib as mpl
from pylab import *
from qutip import *
from matplotlib import cm
import imageio


def animate_bloch(states, name,  duration=0.1, save_all=False):
    b = Bloch()
    b.vector_color = ['k']
    b.view = [-40, 30]
    images = []
    try:
        length = len(states)
    except:
        length = 1
        states = [states]
    ## normalize colors to the length of data ##
    nrm = mpl.colors.Normalize(0, length)
    colors = cm.cool(nrm(range(length)))  # options: cool, summer, winter, autumn etc.

    ## customize sphere properties ##
    b.point_color = list(colors)  # options: 'r', 'g', 'b' etc.
    b.point_marker = ['o']
    b.point_size = [30]

    for i in range(length):
        b.clear()
        b.add_states(states[i])
        b.add_states(states[:(i + 1)], 'point')
        if save_all:
            b.save(dirc='tmp')  # saving images to tmp directory
            filename = "tmp/bloch_%01d.png" % i
        else:
            filename = 'temp_file.png'
            b.save(filename)
        images.append(imageio.imread(filename))
    imageio.mimsave(name, images, duration=duration)

#rotation in x
def rot_x(theta):
    rotation = (np.cos(theta/2) * eye(2) - (1j*  np.sin(theta/2)*sigmax()))
    return rotation

#rotation in y
def rot_y(theta):
    rotation = (np.cos(theta / 2) * eye(2)) - (1j * np.sin(theta / 2) * sigmay())
    return rotation

#rotation in z
def rot_z(theta):
    rotation = np.cos(theta / 2) * eye(2) - 1j * np.sin(theta / 2) * sigmaz()
    return rotation

pz = basis(2,0)
mz = basis(2,1)
px = 1/sqrt(2)*(pz + mz)
mx = 1/sqrt(2)*(pz - mz)
py = 1/sqrt(2)*(pz + 1j*mz)
my = 1/sqrt(2)*(pz - 1j*mz)

print(rot_x(pi) * basis(2,0))

states = []
thetas = linspace(0,pi,21)
for theta in thetas:
    states.append(( (cos(-theta/2)*basis(2,0) - 1j* sin(-theta/2)* sigmax() * basis(2,0))).unit())

animate_bloch(states, 'xrotation.gif', duration=0.1, save_all=False)

states_y = []
thetas = linspace(0,pi,21)
for theta in thetas:
    states_y.append((cos(theta/2)*basis(2,0) - 1j* sin(-theta/2)* sigmay() * basis(2,0)).unit())

animate_bloch(states_y, 'yrotation.gif', duration=0.1, save_all=False)

states_Hadamard = []
thetas = linspace(0,pi,21)
for theta in thetas:
    state = ( (cos(-theta/2)*basis(2,0) - 1j* sin(-theta/2)* sigmax() * basis(2,0))).unit()
    states_Hadamard.append(state)

inter_state = states_Hadamard[-1]
print(inter_state)
thetas2 = linspace(0,pi/2,21)
#length= np.arange(0,len(inter_state),1)
for theta2 in thetas2:
    final_state = (1j * (cos(-theta2/2)*inter_state - 1j* sin(-theta2/2)* sigmay() * inter_state)).unit()
    states_Hadamard.append(final_state)


animate_bloch(states_Hadamard, 'Hadamardrotation.gif', duration=0.1, save_all=False)


states_z = []
thetas = linspace(0,pi,21)
for theta in thetas:
    states_z.append((cos(theta/2) * px - 1j* sin(theta/2) * sigmaz() * px).unit())

animate_bloch(states_z, 'zrotation.gif', duration=0.1, save_all=False)