import matplotlib.pyplot as plt
import cmath
from functions import *

def main():
    ## Fixed Parameters
    m = 1
    E = np.sqrt(8)
    a = 2
    epsilon = 0.01
    # E = 1 # Just adding this for now; it'll be varied
    N = 2000 # Number of discretized momenta

    # Momenta ranges from k_min = 0 to k_max which depends on 3-body CM Frame Energy (E)
    
    # print(M_phib(E, m, a, N, epsilon))
    # print(cmath.sqrt(triangle_function(E, s_b(m, a), m**2)))
    # q_momentum = q(E, m, a)
    # print(k_max(E, m))
    # print(q_momentum)

    print(B(E, m, a, N, epsilon))

if __name__ == '__main__':
    main()