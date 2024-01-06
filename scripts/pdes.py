from typing import Callable

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['lines.linewidth'] = 0.5

def diffusion_1d(
    coefficient: float, 
    left_bc_type: str, right_bc_type: str, 
    left_bc_value: float, right_bc_value: float, 
    initial_condition: Callable, 
    domain_size: float, time_interval: list[float, float], 
    domain_nodes: int, time_steps: int
):

    # DISCRETIZE THE PDE
    # Generate the space domain and split it into Nx nodes.
    x = np.linspace(0, domain_size, domain_nodes)
    # Define the step size in space.
    dx = domain_size / (domain_nodes - 1)

    # Generate the time domain and split it into Nt+1 intervals.
    t = np.linspace(time_interval[0], time_interval[1], time_steps)
    # Define the step size in time.
    dt = (time_interval[1] - time_interval[0]) / (time_steps - 1)

    # Initialize the solution matrix.
    u = np.zeros((domain_nodes, time_steps))
    # Apply the initial condition 
    u[:, 0] = initial_condition(x)

    # Initialize the coefficient matrix.
    D = np.zeros((domain_nodes, domain_nodes))
    C = coefficient * dt / dx**2
    for i in range(1, domain_nodes - 1):
        D[i, i - 1] = C
        D[i, i] = 1 - 2 * C
        D[i, i + 1] = C

    # Apply the boundary conditions.
    # Generate a forcing vector for inhomogeneous Neumann B.C.s.
    F = np.zeros(domain_nodes)
    # Apply the left boundary condition
    match left_bc_type:
        case 'D': 
            u[0, 1:] = left_bc_value
            D[0, 0] = 1
        case 'N':
            F[0] = -2 * left_bc_value * coefficient * dt / dx
            D[0, 0] = 1 - 2 * C
            D[0, 1] = 2 * C
    # Apply the right boundary condition
    match right_bc_type:
        case 'D':
            u[-1, :] = right_bc_value
            D[-1, -1] = 1
        case 'N':
            F[-1] = 2 * right_bc_value * coefficient * dt / dx
            D[-1, -1] = 1 - 2 * C
            D[-1, -2] = 2 * C

    # SOLVE THE PDE
    # Step through each time interval and solve for node values.
    for i in range(1, time_steps):
        u[:, i] = D @ u[:, i - 1] + F
    
    return x, t, u


def wave_1d(
    coefficient: float, 
    left_bc_type: str, right_bc_type: str, 
    left_bc_value: float, right_bc_value: float, 
    first_initial_condition: Callable, 
    second_initial_condition: Callable,
    domain_size: float, time_interval: list[float, float], 
    domain_nodes: int, time_steps: int
):
    # DISCRETIZE THE PDE
    # Generate the space domain and split it into Nx nodes.
    x = np.linspace(0, domain_size, domain_nodes)
    # Define the step size in space.
    dx = domain_size / (domain_nodes - 1)

    # Generate the time domain and split it into Nt+1 intervals.
    t = np.linspace(time_interval[0], time_interval[1], time_steps)
    # Define the step size in time.
    dt = (time_interval[1] - time_interval[0]) / (time_steps - 1)

    # Initialize the solution matrix.
    u = np.zeros((domain_nodes, time_steps))
    # Apply the first initial condition.
    u[:, 0] = first_initial_condition(x)

    # Compute the Courant number.
    r = coefficient * dt / dx

    # Initialize a coefficient matrix for the first time step.
    D0 = np.zeros((domain_nodes, domain_nodes))
    for i in range(1, domain_nodes - 1):
        D0[i, i - 1] = 1/2 * r**2
        D0[i, i] = 1 - r**2
        D0[i, i + 1] = 1/2 * r**2

    # Initialize coefficient matrices for the remaining 3-level time steps.
    D = np.zeros((domain_nodes, domain_nodes))
    E = np.zeros((domain_nodes, domain_nodes))
    for i in range(1, domain_nodes - 1):
        D[i, i - 1] = r**2
        D[i, i] = 2 * (1 - r**2)
        D[i, i + 1] = r**2
        E[i, i] = -1

    # Apply the boundary conditions.
    # Generate a forcing vector for inhomogeneous Neumann B.C.s.
    F = np.zeros(domain_nodes)
    # Apply the left boundary condition
    match left_bc_type:
        case 'D': 
            u[0, 1:] = left_bc_value
            D0[0, 0] = 1
            D[0, 0] = 1
        case 'N':
            F[0] = -2 * left_bc_value * r**2 * dx
            D0[0, 0] = 1 - r**2
            D0[0, 1] = r**2
            D[0, 0] = 2 * (1 - r**2)
            D[0, 1] = 2 * r**2
            E[0, 0] = -1
    # Apply the right boundary condition
    match right_bc_type:
        case 'D':
            u[-1, :] = right_bc_value
            D0[-1, -1] = 1
            D[-1, -1] = 1
        case 'N':
            F[-1] = 2 * right_bc_value * r**2 * dx
            D0[-1, -1] = 1 - r**2
            D0[-1, -2] = r**2
            D[-1, -1] = 2 * (1 - r**2)
            D[-1, -2] = 2 * r**2
            E[-1, -1] = -1
    
    # Implement the derivative initial condition.
    G = dt * second_initial_condition(x)
    u[:, 1] = D0 @ u[:, 0] + F + G

    # SOLVE THE PDE
    # Step through each time interval and solve for node values.
    for t in range(2, time_steps - 1):
        u[:, t] = D @ u[:, t - 1] + E @ u[:, t - 2] + F

    return x, t, u


def laplace_2d(
    left_bc_type: str, right_bc_type: str, 
    top_bc_type: str, bottom_bc_type: str,
    left_bc: Callable, right_bc: Callable, 
    top_bc: Callable, bottom_bc: Callable,
    domain_size_x, domain_size_y,
    domain_nodes_x, domain_nodes_y,
    target_error,
):
    # Initialize the mesh.
    x = np.linspace(0, domain_size_x, domain_nodes_x)
    y = np.linspace(0, domain_size_y, domain_nodes_y)

    # Compute step sizes.
    dx = domain_size_x / (domain_nodes_x - 1)
    dy = domain_size_y / (domain_nodes_y - 1)

    # Add fictitious points.
    x = np.concatenate(([x[0] - dx/2], x, [x[-1] + dx/2]))
    y = np.concatenate(([x[0] - dy/2], y, [y[-1] + dy/2]))

    # Initialize the field on the grid.
    u = np.zeros((domain_nodes_x + 2, domain_nodes_y + 2))

    # Apply the boundary conditions.
    match left_bc_type:
        case 'D':
            u[1, :] = left_bc(x,y)
        case 'N':
            u[0, :] = u[2, :] - 2 * dx * left_bc(x,y)
    match right_bc_type:
        case 'D':
            u[-2, :] = right_bc(x,y)
        case 'N':
            u[-1, :] = u[-3, :] - 2 * dx * right_bc(x,y)
    match top_bc_type:
        case 'D':
            u[:, -2] = top_bc(x,y)
        case 'N':
            u[:, -1] = u[:, -3] - 2 * dy * np.transpose(top_bc(x,y))
    match bottom_bc_type:
        case 'D':
            u[:, 1] = bottom_bc(x,y)
        case 'N':
            u[:, 0] = u[:, 2] - 2 * dy * np.transpose(bottom_bc(x,y))
    
    # Initialize the iteration matrices.
    count = 0
    while True:
        count += 1
        old_u = u.copy()

        for i in range(1, domain_nodes_x + 1):
            for j in range(1, domain_nodes_y + 1):
                # The discrete Laplace equation on a finite difference stencil
                u[i, j] = (
                    (
                        (old_u[i + 1, j] + old_u[i - 1, j]) * dy**2
                        + (old_u[i, j + 1] + old_u[i, j - 1]) * dx**2
                    )
                    / (2 * dx**2 + 2 * dy**2)
                )
        
        # Update the BCs
        match left_bc_type:
            case 'D':
                u[1, :] = left_bc(x,y)
            case 'N':
                u[0, :] = u[2, :] - 2 * dx * left_bc(x,y)
        match right_bc_type:
            case 'D':
                u[-2, :] = right_bc(x,y)
            case 'N':
                u[-1, :] = u[-3, :] - 2 * dx * right_bc(x,y)
        match top_bc_type:
            case 'D':
                u[:, -2] = top_bc(x,y)
            case 'N':
                u[:, -1] = u[:, -3] - 2 * dy * np.transpose(top_bc(x,y))
        match bottom_bc_type:
            case 'D':
                u[:, 1] = bottom_bc(x,y)
            case 'N':
                u[:, 0] = u[:, 2] - 2 * dy * np.transpose(bottom_bc(x,y))
        
        # Compute the norm of the error
        E = np.sum(np.square(u - old_u))

        # Stop iterating if the solution has converged.
        if E <= target_error:
            print(E)
            break
    
    return x[1:-1], y[1:-1], u[1:-1, 1:-1], count


def plot_2d(x: np.ndarray, u: np.ndarray, num_lots = 50):
    _, Nt = u.shape

    for n in range(0, Nt, Nt // num_lots):
        plt.plot(x, u[:, n])

    plt.xlabel('$x$')
    plt.ylabel('$u(x,t)$')
    plt.show()


def plot_3d(
    x: np.ndarray, t: np.ndarray, u: np.ndarray,
    view: tuple[float, float, float]
):
    x, t = np.meshgrid(x, t)

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    surf = ax.plot_surface(x, t, np.transpose(u), cmap=cm.coolwarm, 
                           linewidth=0,
                           )
    ax.view_init(*view)
    plt.show()