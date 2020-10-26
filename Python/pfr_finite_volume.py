# -*- coding: utf-8 -*-
"""
Idealized one-dimensional tubular reactor with first order irreversible 
reaction, convection and diffusion.

@author: Ivo Filot
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    # set system variables
    N = 20         # number of control volumes
    Pe = 1e3    # Peclet number (do not go above 1000!)
    Da = 1.8      # Damköhler number
    
    # build matrix-vector system
    A,b = build_system(N, Pe, Da)
       
    # solve system
    u = np.linalg.solve(A,b)
    
    # plot result
    plt.figure(dpi=150, figsize=(10,4))
    z1 = np.linspace(0.5 / N, 1 - 0.5/N, N, endpoint=True)
    plt.plot(z1, u, 'o', label = 'Numerical integration')
    z2 = np.linspace(0, 1, 100)
    plt.plot(z2, analytical_solution_laplace(z2, Pe, Da), '--', label = 'Analytical solution')
    plt.grid(linestyle='--', color='grey')
    plt.xlabel('Dimensionless axial length $z$ [-]')
    plt.ylabel('Dimensionless concentration $u$ [-]')
    plt.ylim(0,1)
    plt.xlim(0,1)
    plt.legend()
    plt.title('Da = %0.1f, Pe = %0.1f, N=%i' % (Da, Pe, N))
    plt.tight_layout()
    plt.show()

    
def build_system(N, Pe, Da):
    A = np.zeros((N,N)) # build matrix
    b = np.zeros(N)     # build vector
    
    # set some auxiliary variables
    dz = 1.0 / N        # length of a control volume
    beta = 1.0 / (Pe * dz)
    Sp = Da * dz        # source/sink term
    Su = -1.0           # auxiliary boundary value term
    
    # populate center values in matrix
    for i in range(1,N-1):
        A[i,i-1] = beta + 0.5               # aw
        A[i,i+1] = beta - 0.5               # ae
        A[i,i] = -A[i,i-1] - A[i,i+1] - Sp  # ap
    
    # populate boundary values in matrix
    A[0,1] = beta - 0.5                     # ae
    A[0,0] = -A[0,1] - Sp + Su              # ap
    A[N-1, N-2] = beta + 0.5                # aw
    A[N-1, N-1] = -A[N-1, N-2] - Sp         # ap
    
    # populate solution vector
    b[0] = Su
    
    return A,b

def analytical_solution(z, Pe, Da):
    """
    Conventional formulation of the exact solution to the one-dimensional
    convection-diffusion-reaction tubular reactor with first-order
    irreversible kinetics.
    """
    beta = np.sqrt(1 + 4.0 * Da / Pe)
    
    nom = np.exp(0.5 * Pe * (1.0 - beta) * z) + (beta - 1.0) / (beta + 1.0) * np.exp(0.5 * Pe * ((1.0 + beta) * z - 2.0 * beta))
    
    denom = 0.5 * (beta + 1.0) - (beta - 1.0)**2 / (2 * (beta + 1.0)) * np.exp(-Pe * beta)
    
    return nom / denom

def analytical_solution_laplace(z, Pe, Da):
    """
    Alternative formulation of the exact solution to the one-dimensional
    convection-diffusion-reaction tubular reactor with first-order
    irreversible kinetics. This solution is based on the Laplace transformation
    of the problem.
    """
    beta = np.sqrt(1 + 4.0 * Da / Pe)
    
    nom = 2.0 * np.exp(0.5 * Pe * z) * (beta * np.cosh(0.5 * Pe * beta * (1.0 - z)) + np.sinh(0.5 * Pe * beta * (1.0 - z)))
    
    denom = (beta**2 + 1.0) * np.sinh(0.5 * Pe * beta) + 2 * beta * np.cosh(0.5 * Pe * beta)
    
    return nom / denom
    
if __name__ == '__main__':
    main()