import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.fftpack import fft2, ifft2
from scipy.integrate import solve_ivp
from scipy.linalg import kron
from matplotlib.animation import FuncAnimation

def lambda_omega(U, V):
    A22 = U**2 + V**2
    lambda_A = 1 - A22
    omega_A = -beta * A22
    return lambda_A, omega_A

def spc_rhs(t, cond, nx, ny, N, KX, KY, K):
    # Split the state vector
    u_temp = cond[:N] 
    v_temp = cond[N:]  

    # Reshape into 2D arrays
    ut = u_temp.reshape((nx, ny))
    vt = v_temp.reshape((nx, ny))

    # Transform back to real space
    U = np.real(ifft2(ut))
    V = np.real(ifft2(vt))

    # Calculate Lambda and Omega
    l, w = lambda_omega(U, V)

    # Reaction terms
    reaction_U = fft2(l * U - w * V)
    reaction_V = fft2(w * U + l * V)
    
    # Diffusion terms
    diffusion_U = -D1 * K * ut
    diffusion_V = -D2 * K * vt

    # Combine terms and return
    dU_hat = reaction_U + diffusion_U
    dV_hat = reaction_V + diffusion_V

    return np.hstack([dU_hat.flatten(), dV_hat.flatten()])

def cheb(N):
	if N==0: 
		D = 0.; x = 1.
	else:
		n = np.arange(0,N+1)
		x = np.cos(np.pi*n/N).reshape(N+1,1) 
		c = (np.hstack(( [2.], np.ones(N-1), [2.]))*(-1)**n).reshape(N+1,1)
		X = np.tile(x,(1,N+1))
		dX = X - X.T
		D = np.dot(c,1./c.T)/(dX+np.eye(N+1))
		D -= np.diag(np.sum(D.T,axis=0))
	return D, x.reshape(N+1)

def cheb_scaled(n, L):
    """
    Generate scaled Chebyshev differentiation matrix and nodes.
    """
    D, x = cheb(n)  # Chebyshev differentiation matrix and nodes on [-1, 1]
    x = L * x  # Scale nodes to [-L, L]
    D = D / L  # Scale differentiation matrix
    return D, x

def cheb_rhs(t, cond, n, L):
    u = cond[:(n+1)**2].reshape((n+1, n+1))
    v = cond[(n+1)**2:].reshape((n+1, n+1))

    # Reaction terms
    l, w = lambda_omega(u, v)
    reaction_u = l * u - w * v
    reaction_v = w * u + l * v
    
    # Diffusion terms
    diffusion_u = D1 * np.dot(L, u.flatten()).reshape(u.shape)
    diffusion_v = D2 * np.dot(L, v.flatten()).reshape(v.shape)
    
    # Combine and return as 1D array
    du_dt = reaction_u + diffusion_u
    dv_dt = reaction_v + diffusion_v
    return np.hstack([du_dt.flatten(), dv_dt.flatten()])

# Define parameters
tspan = np.arange(0, 10.5, 0.5)
t_eval = tspan  # Times at which to store solution
t_range = (tspan[0], tspan[-1])
beta = 10
Lx, Ly = 20, 20
nx, ny = 64, 64
N = nx * ny
m = 3
D1 = 0.1
D2 = 0.1

# Define spatial domain and initial conditions
x2 = np.linspace(-Lx/2, Lx/2, nx + 1)
x = x2[:nx]
y2 = np.linspace(-Ly/2, Ly/2, ny + 1)
y = y2[:ny]
X, Y = np.meshgrid(x, y)
r = np.sqrt(X**2 + Y**2)
theta = np.angle(X + 1j * Y)
u = np.tanh(r) * np.cos(m * theta - r)
v = np.tanh(r) * np.sin(m * theta - r)

'''
Part a
'''
kx = (2 * np.pi / Lx) * np.concatenate((np.arange(0, nx/2), np.arange(-nx/2, 0)))
kx[0] = 1e-6
ky = (2 * np.pi / Ly) * np.concatenate((np.arange(0, ny/2), np.arange(-ny/2, 0)))
ky[0] = 1e-6
KX, KY = np.meshgrid(kx, ky)
K = KX**2 + KY**2

# Solve the ODE and plot the results
ut0 = fft2(u)
vt0 = fft2(v)
initial = np.concatenate((ut0.flatten(), vt0.flatten()))
sol_spc = solve_ivp(spc_rhs, t_range, initial, t_eval=t_eval, args=(nx, ny, N, KX, KY, K))

'''
Part b
'''
n = 30
L = 20
Lx, Ly = L/2, L/2
D, x = cheb_scaled(n, Lx)
D[0, :] = 0  # No-flux at the left boundary
D[-1, :] = 0  # No-flux at the right boundary

D22= np.dot(D, D)/Lx**2
I = np.eye(len(D22))
L = kron(I, D22) + kron(D22, I)
X, Y = np.meshgrid(x, x)
r = np.sqrt(X**2 + Y**2)
theta = np.angle(X + 1j * Y)
u = np.tanh(r) * np.cos(m * theta - r)
v = np.tanh(r) * np.sin(m * theta - r)

initial_cheb = np.hstack([u.flatten(), v.flatten()])
sol_cheb = solve_ivp(cheb_rhs, t_range, initial_cheb, t_eval=t_eval, args=(n, L))

'''
Plots
'''
# Animation for Fourier Solutions
U_solutions = [np.real(ifft2(sol_spc.y[:N, i].reshape((nx, nx)))) for i in range(len(t_eval))]
V_solutions = [np.real(ifft2(sol_spc.y[N:, i].reshape((nx, nx)))) for i in range(len(t_eval))]

# Create the figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Display the initial frame for both U and V
im_U = ax1.imshow(U_solutions[0], cmap='viridis', vmin=-1, vmax=1)
ax1.set_title("U Field")
ax1.set_xlabel("x")
ax1.set_ylabel("y")

im_V = ax2.imshow(V_solutions[0], cmap='plasma', vmin=-1, vmax=1)
ax2.set_title("V Field")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax_time = fig.add_axes([0.5, 0.92, 0.5, 0.05])  # [left, bottom, width, height]
ax_time.axis('off')  # Turn off the axis for time display

# Time text display outside the animation
time_text = ax_time.text(0.5, 0.5, '', ha='center', va='center', color='black', fontsize=12)

# Animation update function
def update_frame(i):
    im_U.set_array(U_solutions[i])  # Update the U field
    im_V.set_array(V_solutions[i])  # Update the V field
    time_text.set_text(f'Time = {t_eval[i]:.2f}')  # Update the time text outside the animation
    return [im_U, im_V, time_text]

# Create the animation
ani = animation.FuncAnimation(fig, update_frame, frames=len(t_eval), interval=200, blit=True)
# ani.save("unstable-three-spiral.mp4", writer="ffmpeg", fps=15, dpi=500)

# Show the animation
plt.show()

# Animation for Chebychev Solutions
# U_solutions = [sol_cheb.y[:(n+1)**2, i].reshape((n+1, n+1)) for i in range(len(t_eval))]
# V_solutions = [sol_cheb.y[(n+1)**2:, i].reshape((n+1, n+1)) for i in range(len(t_eval))]

# # Create the figure with two subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# # Display the initial frame for both U and V
# im_U = ax1.imshow(U_solutions[0], cmap='viridis', vmin=-1, vmax=1)
# ax1.set_title("U Field")
# ax1.set_xlabel("x")
# ax1.set_ylabel("y")

# im_V = ax2.imshow(V_solutions[0], cmap='plasma', vmin=-1, vmax=1)
# ax2.set_title("V Field")
# ax2.set_xlabel("x")
# ax2.set_ylabel("y")

# # Create the animation
# ani = animation.FuncAnimation(fig, update_frame, frames=len(t_eval), interval=500, blit=True)

# # Show the animation
# plt.show()
