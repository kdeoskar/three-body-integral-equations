#fitting:
#create a list as ydata and N_list as the xdata list
#Define the fitting models
#params, covariance = curve_fit(nonlinear_func, x_data, y_data)
#for linear model, α = params[0] 
#for quadratic model, α, β = params[0], params[1]
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as fitter

#The list of N-data
N_list = [1000 + 500 * i for i in range(11)]   #the elements of N_list can be changed later

#The lists of different y-data
Re_M_data = []
Im_M_data = []

for i in range(11):
    Re_M_data.append(M_phib_bound_state_value(E, m, a, N[i], epsilon).real)
    Im_M_data.append(M_phib_bound_state_value(E, m, a, N[i], epsilon).imag)
    
#Linear model
def linear(x, a0, a1):
    return a0 - a1 / x
#fitting of Re_M_data
par_linear_re, cov_linear_re = curve_fit(linear, N_list, Re_M_data)
a0_re = par_linear_re[0]
a1_re = par_linear_re[1]   #The real part of fitting para α
#fitting of Im_M_data
par_linear_im, cov_linear_im = curve_fit(linear, N_list, Im_M_data)
a0_im = par_linear_im[0]
a1_im = par_linear_im[1]   #The imaginary part of fitting para α

#Quadratic model
def quadratic(x, b0, b1, b2):
    return b0 - b1 / x - b2 / (x**2)
#fitting of Re_M_data
par_quad_re, cov_quad_re = curve_fit(quadratic, N_list, Re_M_data)
b0_re = par_quad_re[0]
b1_re, b2_re = par_quad_re[1], par_quad_re[2]   #The real part of fitting para α, β
#fitting of Im_M_data
par_quad_im, cov_quad_im = curve_fit(quadratic, N_list, Im_M_data)
b0_im = par_quad_im[0]
b1_im, b2_im = par_quad_im[1], par_quad_im[2]   #The imaginary part of fitting para α, β

print("The fitting parameter α of linear model is: {0:6.3f}+{0:6.3f}i".format(a1_re, a1_im))
print("The fitting parameters α, β of quadratic model are: α={0:6.3f}+{0:6.3f}i, β={0:6.3f}+{0:6.3f}i".format(b1_re, b1_im, b2_re, b2_im))

#The M_phib_inverse calculated from the fitting formulae
M_phib_inv_linear_im = []
M_phib_inv_quad_im = []

for i in range(11):
    M_phib_inv_linear_im.append((1 / ((a0_re + 1j*a0_im) - (a1_re + 1j*a1_im) / N[i])).imag)
    M_phib_inv_quad_im.append((1 / ((b0_re + 1j*b0_im) - (b1_re + 1j*b1_im) / N[i] - (b2_re + 1j*b2_im) / (N[i]**2))).imag)
    

#The calculation of Δρ
delta_rho_linear = []
delta_rho_quad = []

for i in range(11):
    delta_rho_linear.append(np.abs(M_phib_inv_linear_im[i] - rho_phib(E, m, a) / rho_phib(E, m, a)) * 100)
    delta_rho_quad.append(np.abs(M_phib_inv_quad_im[i] - rho_phib(E, m, a) / rho_phib(E, m, a)) * 100)
    
#Plotting
