#2D wave equation simulation by Dylan G.

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rcParams["figure.figsize"] = [7.50, 6.50]
plt.rcParams["figure.autolayout"] = True

width = 3
height = 3
resolution = 20

h = 1 / resolution
c = resolution ** 2 / 3000
X = np.arange(-width / 2, width / 2 + h, h)
Y = np.arange(-height / 2, height / 2 + h, h)
X, Y = np.meshgrid(X, Y)

m = width * resolution + 1
n = height * resolution + 1

#Get TST matrices for building block Toeplitz matrices
def get_tsts(z):
    dia = [0] + [1] * (z - 2)
    return sp.spdiags(dia, 0, z, z), sp.spdiags([dia] * 2, [-1,1], z, z).T

#Get forward/backward finite difference coeff's for von Neumann boundary conditions
def vn_fdc(arr, k):
    arr[: k] = (arr[ k :  k*2] * 4 - arr[ k*2 :  k*3]) / 3
    arr[-k:] = (arr[-k*2 : -k] * 4 - arr[-k*3 : -k*2]) / 3
    return arr

#Create block matrices of finite difference coeff's for stepping forward through time
def generate_matrices(vn):
    Im0, Im1 = get_tsts(m)
    A0 = sp.lil_matrix((2 - 4 * c) * Im0 + c * Im1, dtype=float)
    A1 = sp.lil_matrix(c * Im0, dtype=float)
    B0 = sp.lil_matrix(Im0, dtype=float)

    if vn:
        A0 = vn_fdc(A0, 1)
        A1 = vn_fdc(A1, 1)
        B0 = vn_fdc(B0, 1)

    In0, In1 = get_tsts(n)
    A = sp.lil_matrix(sp.kron(In0, A0) + sp.kron(In1, A1))
    B1 = sp.lil_matrix(In0, dtype=float)

    if vn:
        A = vn_fdc(A, m)
        B1 = vn_fdc(B1, 1)

    A = sp.csr_matrix(A)
    B = sp.csr_matrix(sp.kron(B1, B0))
    return A, B

A,B = generate_matrices(True)

##Create matrix sparsity diagram and exit
#plt.spy(A, marker='.', markersize=0.8)
#plt.show()
#exit()

tmax = 200
Z = [np.zeros(n * m)]
Z += [np.zeros(n * m)]
r2 = (X-0.5) ** 2 + (Y-0.3) ** 2
#Z += [0.2 * (np.exp(-20 * r2)).flatten()]# * (1-20*r2)).flatten()]

for t in range(2, tmax):
    Z += [ A * Z[t-1] - B * Z[t-2] ]

    ##Add periodic wave source
    #if vn:
    Z[t] += 0.01*np.exp(-30 * ((X-0.2) ** 2 + (Y-0.5) ** 2)).flatten() * np.sin((t-2)/8)

def change_plot(j, Z, plot):
   plot[0].remove()
   plot[0] = ax.plot_surface(X, Y, Z[j].reshape((n,m)), cmap='viridis', linewidth=0.5, edgecolor = 'k')
   ax.set_zlim(-1, 1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot = [ax.plot_surface(X, Y, Z[0].reshape((n,m)), cmap='viridis', linewidth=0.5, edgecolor = 'k')]
#(plot var is put inside array so that it can be passed into change_plot by reference)

lim = max(width/2, height/2)
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_zlim(-1, 1)

ani = animation.FuncAnimation(fig, change_plot, tmax, fargs=(Z, plot), interval=10)
ax.axis('off')
plt.show()
