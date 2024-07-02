import numpy as np
import matplotlib.pyplot as plt

## Fixed Parameters
m = 1
a = 5

N = 1000 # Number of discretized momenta

# Momenta ranges from k_min = 0 to k_max which depends on 3-body CM Frame Energy (E)
k_min = 0


'''
Helper functions:
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

