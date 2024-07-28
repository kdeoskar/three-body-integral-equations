import numpy as np
import cmath
import math

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
#TODO: the commented out code is actually k_max SQUARED
def k_max(E:float, m:float) -> float:
    #return ((E**2-m**2)/ (2*E))**2
    return (E**2-m**2)/ (2*E)

def k_min():
    return 0

# Function to calculate 
def delta_k(E:float, m:float, N:int) -> float:
    return (k_max(E, m) - k_min()) / N

# Generates list of momenta for energy E and N number of discrete points
# Exclude k_min because it is zero
def momenta(E:float, m:float, N:int):
    momenta_array = np.linspace(delta_k(E, m, N), k_max(E, m), N, endpoint = True) 
    return momenta_array

'''
Helper functions related to Energy:
'''
#Time component of 4-vector k
def omega(m:float, k:float) -> float:
    return np.sqrt(m**2 + k**2)

def alpha(E, m, p, k):
    return (E - omega(m, p) - omega(m, k))**2 - p**2 - k**2 - m**2

# CM Frame energy of particle in the bound state with momentum k
def E_2(E, m, k):
    return (E - omega(m, k))**2 -k**2

# Mandelstam Variable corresponding to 3-momentum k and mass m
def s2k(E, m, k):
    return E_2(E, m, k)**2

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
 
# Removing factors of p and k in denominator for now
def G_S(E, m, p, k, epsilon) -> complex:
    return -H(E, m, p, k) / (4 * p * k ) * cmath.log( (alpha(E, m, p, k) - 2*p*k + 1j*epsilon) /  (alpha(E, m, p, k) + 2*p*k + 1j*epsilon)  ) #1j is the imaginary unit in python syntax

# Bound-state internal scattering amplitude
def M_2(E, m, a, k, epsilon) -> complex:
    return -16*np.pi*cmath.sqrt(s2k(E, m, k) + 1j*epsilon) / (1 / a + 1j*cmath.sqrt(((s2k(E, m, k) + 1j*epsilon)/2)**2 - m**2))

#Calculates the epsilon value according to different eta
def epsilon(E, m, N, η):
    return (η * k_max(E, m)) / (2 * np.pi * N)
    
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
def B_inv(E, m, a, N, epsilon):
    B_inv = np.linalg.inv(B(E, m, a, N, epsilon))
    return B_inv

'''
Function to solve for d_S() (eq 34)
'''
def d_S(E, m, a, p, k, epsilon, N):
    '''Solves for d_S using Brute Force Method'''
    B_inverse = B_inv(E, m, a, N, epsilon)

    # Find the Indices of p and k in the momentum array
    n_p = int(p / delta_k(E, m, N))
    n_k = int(k / delta_k(E, m, N))

    return B_inverse[n_p][n_k] * G_S(E, m, p, k, epsilon)

def d_S_matrix(E, m, a, epsilon, N):
    B_inverse = B_inv(E, m, a, N, epsilon)
    d_matrix = np.zeros((N, N), dtype=complex)
    momenta_array = momenta(E, m, N)
    
    for i in range(len(momenta_array)):
        for j in range(len(momenta_array)):
            k_i = momenta_array[i]
            k_j = momenta_array[j]
            d_matrix[i][j] = - B_inverse[i][j] * G_S(E, m, k_i, k_j, epsilon)

    return d_matrix

###############################################
# The above code gives us d_S(p,k; epsilon, N)
# The next step is to calculate the residue 'g'
###############################################
'''
The pole position in terms of the binding momentum, s_b, (eq 17)
'''
def s_b(m, a):
    return 4*( m**2 - (1 / a)**2 )

'''
Residue g (eq 18)
'''
def residue_g(m, a): 
    return 8 * np.sqrt(2*np.pi * np.sqrt(s_b(m, a)) * (1/a))


###############################################
# Finding the momentum corresponding to the bound state pole s_b, denoted by q
###############################################
'''
Kallen Triangle Function 
'''
def triangle_function(x, y, z):
    return x**2 + y**2 + z**2 - 2*( x*y + y*z + x*z )

'''
"on shell" value of k as q (momentum corresponding bound state pole) (eq 19)
'''
def q(E, m, a) -> complex:
    return (1 / (2*E)) * cmath.sqrt(triangle_function(E**2, s_b(m, a), m**2))

'''
Two body phase space between bound state and spectator (eq 26)
'''
def rho_phib(E, m, a) -> complex:
    return ( q(E, m, a) / (8*np.pi* E ))


################################################
# Now, we just need to take appropriate limits
################################################

'''
Numerical evaluation of the limit to the bound state pole (eq 24)
'''
def M_phib(E, m, a, N, epsilon) -> complex:
    q_momentum = q(E, m, a).real
    return residue_g(m, a)**2 * d_S(E, m, a, q_momentum, q_momentum, epsilon, N)

def M_phib_matrix(E, m, a, N, epsilon):
    return residue_g(m, a)**2 * d_S_matrix(E, m, a, epsilon, N) 

def M_phib_bound_state_value(E, m, a, N, epsilon) -> complex:
    q_momentum = q(E, m, a).real
    q_index = int(q_momentum / delta_k(E, m, N) )
    val = M_phib_matrix(E, m, a, N, epsilon)[q_index][q_index]
    # print("Value of M_phib at the bound state is : ", val)
    # return M_phib_matrix(E, m, a, N, epsilon)[q_index][q_index]
    return val

def M_phib_inv(E, m, a, N, epsilon):
    return np.linalg.inv(M_phib_matrix(E, m, a, N, epsilon)) # THIS IS ERRORING OUT 

def M_phib_inv_bound_state_value(E, m, a, N, epsilon) -> complex:
    q_momentum = q(E, m, a).real
    q_index = int(q_momentum / delta_k(E, m, N) )
    return M_phib_inv(E, m, a, N, epsilon)[q_index][q_index]

################################################
# To calculate uncertainty Delta rho_{\varphi b}
# and re, im parts of rho * M
################################################

def Re_rho_M(E, m, a, N, epsilon):
    return rho_phib(E, m, a) *  (M_phib(E, m, a, N, epsilon)).real

def Re_rho_M_matrix(E, m, a, N, epsilon):
    return (rho_phib(E, m, a) * M_phib_bound_state_value(E, m, a, N, epsilon)).real

def Im_rho_M(E, m, a, N, epsilon):
    return (rho_phib(E, m, a).real * M_phib(E, m, a, N, epsilon)).imag

def Im_rho_M_matrix(E, m, a, N, epsilon):
    return (rho_phib(E, m, a).real *  M_phib_inv_bound_state_value(E, m, a, N, epsilon)).imag

def delta_rho(E, m, a, N, epsilon):
    return np.abs( (Im_rho_M_matrix(E, m, a, N, epsilon) - rho_phib(E, m, a)) / rho_phib(E, m, a) ) * 100


# Returns values for Re(rho_{\varphi b} M_{\varphi b}), Im(rho_{\varphi b} M_{\varphi b}), and Delta rho_{\varphi b}
def return_values(E, m, a, N, epsilon):
    M_bound_state_value = M_phib_bound_state_value(E, m, a, N, epsilon)
    M_inv_bound_state_value =  M_phib_inv_bound_state_value(E, m, a, N, epsilon)
    rho_value = rho_phib(E, m, a)
    list = []

    list[0] = (rho_value * M_bound_state_value).real
    list[1] = (rho_value * M_bound_state_value).imag
    list[2] = np.abs( (M_inv_bound_state_value - rho_value) / rho_value) * 100

    return list


def main():
    pass

if __name__ == '__main__':
    main()

