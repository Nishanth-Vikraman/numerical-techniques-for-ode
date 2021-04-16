"""
Created on Fri Apr  9 16:02:16 2021

@author: nishanth
"""

def f1(x, y, *z):
    return (1 / (x + y))

def f2(x, y, z):
    return (x*z - 4*y)

# RK-4 method
def rk4(x0, y0, xn, n):
    # Calculating step size
    h = (xn-x0)/n
    print('\n--------SOLUTION--------')
    print('-------------------------')
    for i in range(n):
        print('x', i+1, ':', x0, '\ny', i+1, ':', y0)
        print('-------------------------')
        k1 = h * (f1(x0, y0))
        print('K1:', k1)
        k2 = h * (f1((x0+h/2), (y0+k1/2)))
        print('K2:', k2)
        k3 = h * (f1((x0+h/2), (y0+k2/2)))
        print('K3:', k3)
        k4 = h * (f1((x0+h), (y0+k3)))
        print('K4:', k4)
        k = (k1+2*k2+2*k3+k4)/6
        y0 = y0 + k
        x0 = x0 + h
    print('x', i+2, ':', x0, '\ny', i+2, ':', y0)
    return y0
    
def rk4_Simultaneous_ODE(x0, y0, z0, xn, n):
    # Calculating step size
    h = (xn-x0)/n
    print('\n--------SOLUTION--------')
    print('-------------------------')
    for i in range(n):
        print('x', i+1, ':', x0, '\ny', i+1, ':', y0, '\nz', i+1, ':', z0)
        print('-------------------------')
        k1 = h * (f1(x0, y0, z0))
        l1 = h * (f2(x0, y0, z0))
        print('K1:', k1)
        print('L1:', l1, '\n')
        k2 = h * (f1((x0+h/2), (y0+k1/2), (z0+l1/2)))
        l2 = h * (f2((x0+h/2), (y0+k1/2), (z0+l1/2)))
        print('K2:', k2)
        print('L2:', l2, '\n')
        k3 = h * (f1((x0+h/2), (y0+k2/2), (z0+l2/2)))
        l3 = h * (f2((x0+h/2), (y0+k2/2), (z0+l2/2)))
        print('K3:', k3)
        print('L3:', l3, '\n')
        k4 = h * (f1((x0+h), (y0+k3), (z0+l3)))
        l4 = h * (f2((x0+h), (y0+k3), (z0+l3)))
        print('K4:', k4)
        print('L4:', l4, '\n')
        k = (k1+2*k2+2*k3+k4)/6
        l = (l1+2*l2+2*l3+l4)/6
        y0 = y0 + k
        z0 = z0 + l
        x0 = x0 + h
    print('x', i+2, ':', x0, '\ny', i+2, ':', y0, '\nz', i+2, ':', z0)
    
#Milne-Simpson Predictor-Corrector method
#Assuming y1 ,y2 and y3(required for milne-simpson formula) are estimated using Fourth- order Runge kutta method
def Milne_Simpson_PC(x, y, x4):
    h = (x4 - x[0])/4
    '''if x[1] is None:
        x[1] = x[0] + h
        y[1] = rk4(x[0], y[0], x[1], 1)
        print(x[1], y[1])
        x[2] = x[1] + h 
        y[2] = rk4(x[1], y[1], x[2], 1)
        x[3] = x[2] + h
        y[3] = rk4(x[2], y[2], x[3], 1)'''
    #Milne Predictor formula
    yp4 = y[0] + 4*h*(2*f1(x[1],y[1]) - f1(x[2],y[2]) + 2*f1(x[3],y[3])) / 3
    fp4 = f1(x4,yp4) #
    print('yp4 = ', y[0], '+ 4 *', h, '* (2 *', f1(x[1],y[1]), '-', f1(x[2],y[2]), '+ 2 *', f1(x[3],y[3]), ') / 3')
    print('    = ',yp4)
    print('fp4 = ',fp4)
    #Simpson Corrector formula
    iterationNumber = 1
    intermediate = y[2] + h*(f1(x[2],y[2]) + 4*f1(x[3],y[3])) / 3
    yc4 = intermediate + h*fp4/3 
    f4 = f1(x4,yc4)
    print('\n---------------------------------------------------------\n')
    print('Iteration Number:', iterationNumber)
    print('yc4 = ', y[0], '+', h, '* (', f1(x[2],y[2]), '+ 4 *', f1(x[3],y[3]), '+', fp4, ') / 3')
    print('    = ', intermediate, '+', h, '*', fp4, '/ 3')
    print('    = ',yc4)
    print('f4 = ',f4)
    while(abs(f4 - fp4) > 0.01):
        fp4 = f4
        iterationNumber += 1
        yc4 = intermediate + h*fp4/3 
        f4 = f1(x4,yc4)
        print('\n---------------------------------------------------------\n')
        print('Iteration Number:', iterationNumber)
        print('yc4 = ', intermediate, '+', h, '*', fp4, '/ 3')
        print('    = ',yc4)
        print('f4 = ',f4)
        
def Adam_PC(x, y, x4):
    h = (x4 - x[0])/4
    '''if x[1] is None:
        x[1] = x[0] + h
        y[1] = rk4(x[0], y[0], x[1], 1)
        print(x[1], y[1])
        x[2] = x[1] + h 
        y[2] = rk4(x[1], y[1], x[2], 1)
        x[3] = x[2] + h
        y[3] = rk4(x[2], y[2], x[3], 1)'''
    #Predictor formula
    yp4 = y[3] + h*(37*f1(x[1],y[1]) - 59*f1(x[2],y[2]) + 55*f1(x[3],y[3]) - 9*f1(x[0], y[0])) / 24
    fp4 = f1(x4,yp4) #
    print('yp4 = ', y[0], '+ 4 *', h, '* (2 *', f1(x[1],y[1]), '-', f1(x[2],y[2]), '+ 2 *', f1(x[3],y[3]), ') / 3')
    print('    = ',yp4)
    print('fp4 = ',fp4)
    #Corrector formula
    iterationNumber = 1
    intermediate = y[3] + h*(-5*f1(x[2],y[2]) + 19*f1(x[3],y[3]) + f1(x[1],y[1])) / 24
    yc4 =  intermediate + h*(9*fp4) / 24
    f4 = f1(x4,yc4)
    print('\n---------------------------------------------------------\n')
    print('Iteration Number:', iterationNumber)
    print('yc4 = ', y[3], '+', h, '* ( 19 *', f1(x[3],y[3]), '- 5 *', f1(x[2],y[2]), '+', f1(x[1],y[1]), '+ 9 *', fp4, ') / 3')
    print('    = ', intermediate, '+', h, '* ( 9 * ', fp4, ') / 24')
    print('    = ',yc4)
    print('f4 = ',f4)
    while(abs(f4 - fp4) > 0.001):
        fp4 = f4
        iterationNumber += 1
        yc4 = intermediate + h*(9*fp4) / 24
        f4 = f1(x4,yc4)
        print('\n---------------------------------------------------------\n')
        print('Iteration Number:', iterationNumber)
        print('yc4 = ', intermediate, '+', h, '* ( 9 * ', fp4, ') / 24')
        print('    = ',yc4)
        print('f4 = ',f4)


'''
# Inputs
print('Enter initial conditions:')
x0 = float(input('x0 = '))
y0 = float(input('y0 = '))

print('Enter calculation point: ')
xn = float(input('xn = '))

print('Enter number of steps:')
step = int(input('Number of steps = '))
'''

#rk4(0,0,1,2)
#rk4_Simultaneous_ODE(0, 3, 0, 0.1, 1)
x = [0.0, 0.2, 0.4, 0.6]
y = [2.0, 2.0933, 2.1755, 2.2532]
Milne_Simpson_PC(x, y, 0.8)
#Adam_PC(x, y, 2)