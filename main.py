import numpy as np
import matplotlib.pyplot as plt
import cmath

## Fixed Parameters
m = 1
a = 5
E = 1 # Just adding this for now; it'll be varied
N = 1000 # Number of discretized momenta

# Momenta ranges from k_min = 0 to k_max which depends on 3-body CM Frame Energy (E)
k_min = 0


'''
Helper functions related to Momentum:
'''

# Function to calculate k_max corresponding to energy E
def k_max(E):
    return ((E**2-m**2)/ (2*E))**2

# Function to calculate 
def delta_k(E, N):
    return (k_max(E) - k_min) / N

# Generates list of momenta for energy E and N number of discrette points
def momenta(E, N):
    momenta_array = np.linspace(k_min, k_max(E), delta_k(E, N), endpoint = True) 
    return momenta_array


'''
Helper functions related to Energy:
'''
#Time component of 4-vector k
def omega(k):
    return np.sqrt(m**2 + k**2)

def alpha(E, p, k):
    return (E - omega(p) - omega(k))**2 - p**2 - k**2 - m**2

# CM Frame energy of particle in the bound state with momentum k
def E_2(E, k):
    return (E - omega(k))**2 -k**2

# Defining cut-off function(s)
def J(x):
    if (x <= 0 ):
        return 0
    elif (x > 1):
        return 1
    else:
        return np.exp(-1/x * np.exp(-1/(1-x)))

def H(p, k):
    return J(E_2(E, p)**2 / (4 * m**2)) * J(E_2(E, k)**2 / (4 * m**2))

def G_S(p, k, epsilon):
    return -H(p,k) / (4 * p * k) * cmath.log( (alpha(E, p,k) - 2*p*k + 1j*epsilon) /  (alpha(E, p,k) + 2*p*k + 1j*epsilon)  ) #1j is the imaginary unit in python syntax

# Bound-state internal scattering amplitude
def M_2(E, k, epsilon):
    return -16*np.pi*(m*a)(E_2(E, k) / m) / (1 + 1j*cmath.sqrt((E_2(E, k)/m)**2 - 1))

# Finds the B matrix and returns its inverse
def B(E, N, epsilon):
    B = []
    momenta_array = momenta(E, N)
    for i in range(N):
        for j in range(N):
            k_i = momenta_array[i]
            k_j = momenta_array[j]
            B[i, j] = np.kron(i, j) + (delta_k(E, N) * k_j**2) / ((2 * np.pi)**2 * omega(k_j)) * G_S(k_i, k_j, epsilon) * M_2(E, k_j, epsilon)

    return B

# Find inverse
def B_inv(E, N, epsilon):
    B_inv = np.linalg.inv(B)
    return B_inv


'''
Function to solve for d_S()
'''
def d_S();