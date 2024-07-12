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
    N = 1000 # Number of discretized momenta

    # Momenta ranges from k_min = 0 to k_max which depends on 3-body CM Frame Energy (E)
    
    # print(M_phib(E, m, a, N, epsilon))
    # print(cmath.sqrt(triangle_function(E, s_b(m, a), m**2)))
    q_momentum = q(E, m, a)
    # print(k_max(E, m))
    # print(bool(k_max(E, m) > q_momentum))
    print(q_momentum)

    start_time = time.time()

    # print(B(E, m, a, N, epsilon))
    print(d_S(E, m, a, q_momentum, q_momentum, epsilon, N))

    print()
    print("Process finished --- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()