import numpy as np
import cmath

# ## Fixed Parameters
# m = 1
# a = 5
# # E = 1 # Just adding this for now; it'll be varied
# N = 1000 # Number of discretized momenta

# # Momenta ranges from k_min = 0 to k_max which depends on 3-body CM Frame Energy (E)
# k_min = 0


'''
Helper functions related to Momentum:
'''

# Function to calculate k_max corresponding to energy E
def k_max(E:int, m:int) -> int:
    return ((E**2-m**2)/ (2*E))**2

def k_min():
    return 0

# Function to calculate 
def delta_k(E:int, m:int, N:int) -> int:
    return (k_max(E, m) - k_min()) / N

# Generates list of momenta for energy E and N number of discrete points
# Exclude k_min because it is zero
def momenta(E:int, m:int, N:int) -> int:
    momenta_array = np.linspace(delta_k(E, m, N), k_max(E, m), N, endpoint = True) 
    return momenta_array


'''
Helper functions related to Energy:
'''
#Time component of 4-vector k
def omega(m:int, k:int) -> int:
    return np.sqrt(m**2 + k**2)

def alpha(E, m, p, k):
    return (E - omega(m, p) - omega(m, k))**2 - p**2 - k**2 - m**2

# CM Frame energy of particle in the bound state with momentum k
def E_2(E, m, k):
    return (E - omega(m, k))**2 -k**2

# Defining cut-off function(s)
def J(x):
    if (x <= 0 ):
        return 0
    elif (x > 1):
        return 1
    else:
        return np.exp(-1/x * np.exp(-1/(1-x)))

def H(E, m, p, k):
    return J(E_2(E, m, p)**2 / (4 * m**2)) * J(E_2(E, m, k)**2 / (4 * m**2))

def G_S(E, m, p, k, epsilon) -> complex:
    return -H(E, m, p, k) / (4 * p * k) * cmath.log( (alpha(E, m, p, k) - 2*p*k + 1j*epsilon) /  (alpha(E, m, p, k) + 2*p*k + 1j*epsilon)  ) #1j is the imaginary unit in python syntax

# Bound-state internal scattering amplitude
def M_2(E, m, a, k, epsilon) -> complex:
    return -16*np.pi*(m*a)*((E_2(E, m, k) + 1j*epsilon) / m) / (1 + 1j*cmath.sqrt(((E_2(E, m, k) + 1j*epsilon)/m)**2 - 1))

# Calculates the B matrix and returns its inverse
def B(E, m, a, N, epsilon):
    B = np.zeros((N, N), dtype=complex)
    momenta_array = momenta(E, m, N)
    for i in range(N):
        for j in range(N):
            k_i = momenta_array[i]
            k_j = momenta_array[j]
            B[i][j] = np.kron(i, j) + (delta_k(E, m, N) * k_j**2) / ((2 * np.pi)**2 * omega(m, k_j)) * G_S(E, m, k_i, k_j, epsilon) * M_2(E, m, a, k_j, epsilon)

    return B

# Find inverse
def B_inv(E, m, N, epsilon):
    B_inv = np.linalg.inv(B(E, m, N, epsilon))
    return B_inv


'''
Function to solve for d_S()
'''
def d_S(E, m, p, k, epsilon, N):
    '''Solves for d_S using Brute Force Method'''
    B_inverse = B_inv(E, m, N, epsilon)

    # Find the Indices of p and k in the momentum array
    n_p = p / delta_k(E, m, N) 
    n_k = k / delta_k(E, m, N)

    return B_inverse[n_p][n_k] * G_S(E, m, p, k, epsilon)

###############################################
# The above code gives us d_S(p,k; epsilon, N)
# The next step is to calculate the residue 'g'
###############################################

def s_b(m, a):
    return 4*( m**2 - (1 / a)**2 )

def residue_g(m, a):
    return 8 * np.sqrt(2*np.pi * np.sqrt(s_b(m, a)) * (1/a))

###############################################
# Finding the momentum corresponding to the bound state pole s_b, denoted by q
###############################################

def triangle_function(x, y, z):
    return x**2 + y**2 + z**2 - 2*( x*y + y*z + x*z )

def q(E, m, a) -> complex:
    return (1 / 2*E) * cmath.sqrt(triangle_function(E**2, s_b(m, a), m**2))

def rho_phib(E, m, a):
    return ( q(E, m, a) / (8*np.pi* E ))

################################################
# Now, we just need to take appropriate limits
################################################

def M_phib(E, m, a, N, epsilon):
    # q = q(E, m, a)
    return residue_g(m, a)**2 * d_S(E, m, q(E, m, a), q(E, m, a), epsilon, N)

def Im_M(E, m, a, N, epsilon):
    pass

def main():
    pass

if __name__ == '__main__':
    main()

