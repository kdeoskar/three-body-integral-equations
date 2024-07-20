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
    N = 2000 # Number of discretized momenta

    # Momenta ranges from k_min = 0 to k_max which depends on 3-body CM Frame Energy (E)
    
    # q_momentum = q(E, m, a).real

    start_time = time.time()

    # print()
    # print(M_phib_matrix(E, m, a, N, epsilon))
    # print()

    # print(Re_rho_M_matrix(E, m, a, N, epsilon))
    # print(Im_rho_M_matrix(E, m, a, N, epsilon))
    # print(delta_rho(E, m, a, N, epsilon))

    # quantities = return_values(E, m, a, N, epsilon)
    # print("Re(rho * M) = ", quantities[0])
    # print("Im(rho * M) = ", quantities[1])
    # print("Delta rho = ", quantities[2])

    '''Random test to see if np.linalg.inv works with complex matrices'''
    # complex_matrix = np.array([
    #     [1e-12+1e-12j, 2e-12+0e-23j],
    #     [3e-12+1e-12j, 2e-12+4e-23j]
    # ], dtype=complex)

    # print(complex_matrix)
    # print(np.linalg.det(complex_matrix))
    # print(complex_matrix[0][0]*complex_matrix[1][1] - complex_matrix[0][1]*complex_matrix[1][0])

    '''Let's see if reducing the scale of the entries makes determinant non-zero'''
    # M_matrix = M_phib_matrix(E,m,a,N,epsilon)
    # print(M_matrix)
    # print(np.linalg.det(M_matrix))

    # M_reduced = (1e-15)*M_matrix
    # print(M_reduced)
    # print(np.linalg.det(M_reduced))

    # M_reduced_inv = np.linalg.inv(M_reduced)
    # print(M_reduced_inv)

    # M_matrix_inv = M_phib_inv(E,m,a,N,epsilon)
    # print(M_matrix_inv)

    # print(d_S_matrix(E, m, a, epsilon, N))

    momenta_array = momenta(E, m, N)
    for p in momenta_array:
        for k in momenta_array:

    #         print(G_S(E, m, p, k, epsilon))
    # print(B(E, m, a, N, epsilon))
    # print(B_inv(E, m, a, N, epsilon))

    # print(rho_phib(E, m, a) * M_phib(E, m, a, N, epsilon))

    print()
    print("Process finished --- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()

