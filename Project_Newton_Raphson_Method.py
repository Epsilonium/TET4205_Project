import numpy as np
import cmath
import math
import pandas as pd
import sympy as sym
"""The Code works, but it keep spinning between values and not converging"""

filename = "Matrix.xlsx"
Z_real = pd.read_excel(filename, sheet_name = 'Z_real', index_col = [0])
Z_imag = pd.read_excel(filename, sheet_name = 'Z_imag', index_col = [0])
Ybus_real = pd.read_excel(filename, sheet_name = 'Y_real', index_col = [0])
Ybus_imag = pd.read_excel(filename, sheet_name = 'Y_imag', index_col = [0])
Ybus = pd.read_excel(filename, sheet_name = 'Ybus', index_col = [0])
Angles = pd.read_excel(filename, sheet_name = 'Angles', index_col = [0])

Ybus_array = np.zeros((5,5), dtype=np.complex_)
Angles_array = np.zeros((5,5))
for x in range(0, len(Ybus_real)):
    for y in range(0, len(Ybus_real)):
        Ybus_array[x][y] = complex(Ybus_real.iat[x,y], Ybus_imag.iat[x, y])
        Angles_array[x][y] = Angles.iat[x, y]


d1, d2, d3, d4, d5, v1, v2, v3, v4, v5, \
p1, p2, p3, p4, p5, q1, q2, q3, q4, q5,\
y11, y12, y13, y14, y15, y21, y22, y23, y24, y25, \
y31, y32, y33, y34, y35, y41, y42, y43, y44, y45, y54, y55, \
th11, th12, th13, th14, th15, th21, th22, th23, th24, th25,\
th31, th32, th33, th34, th35, th41, th42, th43, th44, th45, th54, th55 = sym.symbols('d1, d2, d3, d4, d5, v1, v2, v3, v4, v5, p1, p2, p3, p4, p5, q1, q2, q3, q4, q5, y11, y12, y13, y14, y15, y21, y22, y23, y24, y25, y31, y32, y33, y34, y35, y41, y42, y43, y44, y45, y54, y55, th11, th12, th13, th14, th15, th21, th22, th23, th24, th25, th31, th32, th33, th34, th35, th41, th42, th43, th44, th45, th54, th55', real=True)


uknown_matrix = [d1, d2, d4, d5, v2, v4, v5]
uknown_matrix_val = [0, 0, 0, 0, 1, 1, 1]
Volt_0 = [1, 1, 1, 1, 1]
delta_0 = [0, 0, 0, 0, 0]



P1_spec = 1
P2_spec = 0.6
P3_spec = 0
P4_spec = 0.6
P5_spec = 0.5
Q1_spec = 0
Q2_spec = 0.3
Q3_spec = 0
Q4_spec = 0.2
Q5_spec = 0.4
Known_spec = [P1_spec, P2_spec, P4_spec, P5_spec, Q2_spec, Q4_spec, Q5_spec]
Power = [P1_spec, P2_spec, P3_spec, P4_spec, P5_spec]
Reactive = [Q1_spec, Q2_spec, Q3_spec, Q4_spec, Q5_spec]


P1_0 = v1*v1*y11*sym.cos(th11)+v1*v2*y12*sym.cos(th12-d1+d2)+v1*v4*y14*sym.cos(th14-d1+d4)
P2_0 = v2*v1*y21*sym.cos(th21-d2+d1)+v2*v2*y22*sym.cos(th22)+v2*v3*y23*sym.cos(th23-d2+d3)+v2*v4*y24*sym.cos(th24-d2+d4)
P4_0 = v4*v1*y41*sym.cos(th41-d4+d1)+v4*v2*y42*sym.cos(th42-d4+d2)+v4*v4*y44*sym.cos(th44)+v4*v5*y45*sym.cos(th45-d4+d5)
P5_0 = v5*v4*y54*sym.cos(th54-d5+d4)+v5*v5*y55*sym.cos(th55)
Q2_0 = v2*v1*y21*sym.sin(d2-d1-th21)+v2*v2*y22*sym.sin(-th22)+v2*v3*y23*sym.sin(d2-d3-th23)+v2*v4*y24*sym.sin(d2-d4-th24)
Q4_0 = v4*v1*y41*sym.sin(d4-d1-th41)+v4*v2*y42*sym.sin(d4-d2-th42)+v4*v4*y44*sym.sin(-th44)+v4*v5*y45*sym.sin(d4-d5-th45)
Q5_0 = v5*v4*y54*sym.sin(d5-d4-th54)+v5*v5*y55*sym.sin(-th55)

P_matrix = sym.Matrix([P1_0, P2_0, P4_0, P5_0, Q2_0, Q4_0, Q5_0])
Jacobi = P_matrix.jacobian(uknown_matrix)

def substitute(any_function, volt, bus, theta, delta):
    """Denne funksjonen skal hente inn et uttrykk og bytte ut
    de symbolske verdiene med tallverdier. Først lages en dictionary
    hvor hvert symbol får en verdi. Derreter skal uttrykket byttes basert
    på hvilke symbolske verdier uttrykket har"""
    dictionary = {}
    for x in range(1, 6):
        v = "v" + str(x)
        d = "d" + str(x)
        dictionary[v] = abs(volt[x-1])
        dictionary[d] = delta[x-1]
        for y in range(1, 6):
            ybus = "y" + str(x) + str(y)
            angle = "th" + str(x) + str(y)
            dictionary[ybus] = abs(bus[x-1][y-1])
            dictionary[angle] = theta[x-1][y-1]
    symbols_used = list(any_function.free_symbols)
    for x in range(0, len(symbols_used)):
        any_function = any_function.subs(symbols_used[x], dictionary[str(symbols_used[x])])

    return any_function
"""
P1_calc = substitute(P1_0, Volt_0, Ybus_array, Angles_array, delta_0)
P2_calc = substitute(P2_0, Volt_0, Ybus_array, Angles_array, delta_0)
P4_calc = substitute(P4_0, Volt_0, Ybus_array, Angles_array, delta_0)
P5_calc = substitute(P5_0, Volt_0, Ybus_array, Angles_array, delta_0)
Q2_calc = substitute(Q2_0, Volt_0, Ybus_array, Angles_array, delta_0)
Q4_calc = substitute(Q4_0, Volt_0, Ybus_array, Angles_array, delta_0)
Q5_calc = substitute(Q5_0, Volt_0, Ybus_array, Angles_array, delta_0)
"""

def next_angles(Jacobi_val, Pmatrix_val, delta, specified):
    """This function finds the next values for uknown variables.
    It needs the jacobi matrix with numbers, the known variables with numbers
    the previous values for the uknown variables and known power values.
    """
    next_uknown = np.zeros((7,1))  # create 7X1 array for the next values
    Delta_known = np.zeros((7,1))  # create 7X1 array for P_spec - P_calc
    print(Jacobi_val)
    Jacobi_val = np.array(Jacobi_val).astype('float')  #Convert the Jacobi matrix to an array. Maybe trouble when converting a sympy matrix with numpy
    for i in range(0, len(specified)):
        Delta_known[i] = (specified[i]-Pmatrix_val[i])  #Fill in the array for P_spec - P_calc
    print("Delta P og Q kjente: \n", Delta_known)
    inverse_jacob = np.linalg.inv(Jacobi_val)    #Take the inverse of Jacobian
    Big_delta = inverse_jacob @ Delta_known    #Calculate Delta values of the uknows by multiplying the inverse of jacib with P_spec - P_calc
    for i in range(0, len(delta)):
        print("før ukjente: ", delta[i])
        next_uknown[i] = delta[i] + Big_delta[i]   #Calculate the next values of the uknowns
    print("de neste ukjente:\n", next_uknown)
    return next_uknown

def Newton_rhap(known_matrix, Jacobi_matrix):
    tol = 0.0001
    k = 1
    iteration_P = substitute(known_matrix, Volt_0, Ybus_array, Angles_array, delta_0)
    iteration_J = substitute(Jacobi_matrix, Volt_0, Ybus_array, Angles_array, delta_0)
    updated_volt = [1 for x in range(0, len(Volt_0))]
    updated_delta = [0 for x in range(0, len(delta_0))]
    updated_uknows = uknown_matrix_val.copy()
    iteration = next_angles(iteration_J, iteration_P, updated_uknows, Known_spec)
    while k<10:
        print("Number of iterations: ", k)
        updated_volt = [1, float(iteration[4]), 1, float(iteration[5]), float(iteration[6])]
        updated_delta = [float(iteration[0]), float(iteration[1]), 0.0, float(iteration[2]), float(iteration[3])]
        updated_uknows = iteration.copy()
        print("Oppdatert: ", updated_volt)
        print("før: ", Volt_0)
        iteration_P = substitute(known_matrix, updated_volt, Ybus_array, Angles_array, updated_delta)
        iteration_J = substitute(Jacobi_matrix, updated_volt, Ybus_array, Angles_array, updated_delta)
        iteration = next_angles(iteration_J, iteration_P, updated_uknows, Known_spec)
        k +=1
    return iteration

print(Newton_rhap(P_matrix, Jacobi))
