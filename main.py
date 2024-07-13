import numpy
import matplotlib.pyplot as plt
import time 
import cmath
from functions import *

def main():
    ## Fixed Parameters
    m = 1
    E = np.sqrt(8)
    a = 2
    epsilon = 0.01
    # E = 1 # Just adding this for now; it'll be varied
    N = 95 # Number of discretized momenta

    # Momenta ranges from k_min = 0 to k_max which depends on 3-body CM Frame Energy (E)
    
    q_momentum = q(E, m, a).real

    start_time = time.time()

    print(Re_rho_M_matrix(E, m, a, N, epsilon))
    print(Im_rho_M_matrix(E, m, a, N, epsilon))
    print(delta_rho(E, m, a, N, epsilon))

    print()
    print("Process finished --- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()