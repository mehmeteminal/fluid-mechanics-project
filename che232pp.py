import math
def Re_function(D,rho,u_mean,mu):
    return rho*u_mean*D/mu
def fF_L(Re): #laminar flow friction factor
    return 16/Re
def fF_T(eps,D,Re): #turbulent flow friction factor
    a = 0.269*eps/D
    b = 2.185/Re
    c = 14.5/Re
    return (-1.737*math.log(a-(b*math.log(a+c,math.e)),math.e))**(-2)
def error_function(old_value,new_value):
    return abs(old_value-new_value)*100/old_value
def case1_SI(Q,D,rho,mu,eps,L,deltaz):
    g = 9.81 #m/s^2
    u_mean = 4*Q/(math.pi*D**2)
    Re = Re_function(D,rho,u_mean,mu)
    if Re<=2000:
        fF = fF_L(Re)
    elif Re>4000:
        fF = fF_T(eps,D,Re)
    else:
        print('Warning: Transition Region!!!')
    deltaP = (2*fF*rho*(u_mean**2)*L/D)+ rho*g*deltaz
    return deltaP
def case1_FPS(Q,D,rho,mu,eps,L,deltaz):
    g = 32.17 #ft/s^2
    g_c = 32.2 #(lbm*ft)/(lbf*s^2)
    u_mean = 4*Q/(math.pi*D**2)
    Re = Re_function(D,rho,u_mean,mu)
    if Re<=2000:
        fF = fF_L(Re)
    elif Re>4000:
        fF = fF_T(eps,D,Re)
    else:
        print('Warning: Transition Region!!!')
    deltaP = (2*fF*rho*(u_mean**2)*L/D)+ rho*g*deltaz
    return deltaP/(g_c*144) # 1ft=12in
def case2_SI(deltaP,D,rho,mu,eps,L,deltaz):
    g = 9.81 #m/s^2
    old_Re = 100000
    if old_Re <= 2000:
        fF = fF_L(old_Re)
    elif old_Re > 4000:
        fF = fF_T(eps, D, old_Re)
    else:
        print('Warning: Transition Region!!!')
    u_mean = (D*(deltaP-rho*g*deltaz)/(2*fF*rho*L))**(1/2)
    new_Re = Re_function(D, rho, u_mean, mu)
    tol = 0.001
    while error_function(old_Re,new_Re)>tol:
        old_Re = new_Re
        if old_Re <= 2000:
            fF = fF_L(old_Re)
        elif old_Re > 4000:
            fF = fF_T(eps, D, old_Re)
        else:
            print('Warning: Transition Region!!!')
        u_mean = (D * (deltaP - rho * g * deltaz) / (2 * fF * rho * L)) ** (1 / 2)
        new_Re = Re_function(D, rho, u_mean, mu)
    Q = u_mean*math.pi*(D**2)/4
    return Q
def case2_FPS(deltaP,D,rho,mu,eps,L,deltaz):
    g = 32.17 #ft/s^2
    g_c = 32.2 #(lbm*ft)/(lbf*s^2)
    old_Re = 100000
    if old_Re <= 2000:
        fF = fF_L(old_Re)
    elif old_Re > 4000:
        fF = fF_T(eps, D, old_Re)
    else:
        print('Warning: Transition Region!!!')
    u_mean = (D*(deltaP*g_c*144-rho*g*deltaz)/(2*fF*rho*L))**(1/2) # 1ft=12in
    new_Re = Re_function(D, rho, u_mean, mu)
    tol = 0.001
    while error_function(old_Re,new_Re)>tol:
        old_Re = new_Re
        if old_Re <= 2000:
            fF = fF_L(old_Re)
        elif old_Re > 4000:
            fF = fF_T(eps, D, old_Re)
        else:
            print('Warning: Transition Region!!!')
        u_mean = (D * (deltaP*g_c*144 - rho * g * deltaz) / (2 * fF * rho * L)) ** (1 / 2) # 1ft=12in
        new_Re = Re_function(D, rho, u_mean, mu)
    Q = u_mean*math.pi*(D**2)/4
    return Q
def case3_SI(deltaP,Q,rho,mu,eps,L,deltaz):
    g = 9.81 #m/s^2
    old_D = 10
    u_mean = 4*Q/(math.pi*old_D**2)
    old_Re = Re_function(old_D,rho,u_mean,mu)
    if old_Re <= 2000:
        fF = fF_L(old_Re)
    elif old_Re > 4000:
        fF = fF_T(eps, old_D, old_Re)
    else:
        print('Warning: Transition Region!!!')
    new_D = ((32*fF*rho*(Q**2)*L)/((math.pi**2)*(deltaP-rho*g*deltaz)))**(1/5)
    u_mean = 4*Q/(math.pi*new_D**2)
    new_Re = Re_function(new_D,rho,u_mean,mu)
    tol = 0.001
    while error_function(old_Re,new_Re)>tol:
        old_Re = new_Re
        if new_Re <= 2000:
            fF = fF_L(new_Re)
        elif new_Re > 4000:
            fF = fF_T(eps, new_D, new_Re)
        else:
            print('Warning: Transition Region!!!')
        new_D = ((32 * fF * rho * (Q ** 2) * L) / ((math.pi ** 2)*(deltaP - rho * g * deltaz))) ** (1 / 5)
        u_mean = 4 * Q / (math.pi * new_D ** 2)
        new_Re = Re_function(new_D, rho, u_mean, mu)
    D = new_D
    return D
def case3_FPS(deltaP,Q,rho,mu,eps,L,deltaz):
    g = 32.17 #ft/s^2
    g_c = 32.2 #(lbm*ft)/(lbf*s^2)
    old_D = 10
    u_mean = 4*Q/(math.pi*old_D**2)
    old_Re = Re_function(old_D,rho,u_mean,mu)
    if old_Re <= 2000:
        fF = fF_L(old_Re)
    elif old_Re > 4000:
        fF = fF_T(eps, old_D, old_Re)
    else:
        print('Warning: Transition Region!!!')
    new_D = ((32*fF*rho*(Q**2)*L)/((math.pi**2)*(deltaP*g_c*144-rho*g*deltaz)))**(1/5) # 1ft=12in
    u_mean = 4*Q/(math.pi*new_D**2)
    new_Re = Re_function(new_D,rho,u_mean,mu)
    tol = 0.001
    while error_function(old_Re,new_Re)>tol:
        old_Re = new_Re
        if new_Re <= 2000:
            fF = fF_L(new_Re)
        elif new_Re > 4000:
            fF = fF_T(eps, new_D, new_Re)
        else:
            print('Warning: Transition Region!!!')
        new_D = ((32 * fF * rho * (Q ** 2) * L) / ((math.pi ** 2)*(deltaP*g_c*144 - rho * g * deltaz))) ** (1 / 5) #1ft=12in
        u_mean = 4 * Q / (math.pi * new_D ** 2)
        new_Re = Re_function(new_D, rho, u_mean, mu)
    D = new_D
    return D
case = int(input('Please specify the case 1 2 or 3: '))
unit = input('Please specify the unit SI or FPS: ')
if case ==1:
    if unit == 'SI':
        Q = float(input('Please specify the flow rate (m^3/s): '))
        D = float(input('Please specify the pipe diameter (m): '))
        rho = float(input('Please specify the density of the fluid (kg/m^3): '))
        mu = float(input('Please specify the viscosity of the fluid (Pa*s): '))
        eps = float(input('Please specify the pipe surface roughness (m): '))
        L = float(input('Please specify the pipe length (m): '))
        deltaz = float(input('Please specify the elevation difference (m): '))
        deltap = case1_SI(Q,D,rho,mu,eps,L,deltaz)
        print('(-)delta pressure drop (Pa): ',deltap)
    else :
        Q = float(input('Please specify the flow rate (ft^3/s): '))
        D = float(input('Please specify the pipe diameter (ft): '))
        rho = float(input('Please specify the density of the fluid (lbm/ft^3): '))
        mu = float(input('Please specify the viscosity of the fluid (lbm/(ft s)): '))
        eps = float(input('Please specify the pipe surface roughness (ft): '))
        L = float(input('Please specify the pipe length (ft): '))
        deltaz = float(input('Please specify the elevation difference (ft): '))
        deltap = case1_FPS(Q, D, rho, mu, eps, L, deltaz)
        print('(-)delta pressure drop (psi): ', deltap)
elif case == 2:
    if unit == 'SI':
        deltaP = float(input('Please specify the (-) delta pressure drop (Pa): '))
        D = float(input('Please specify the pipe diameter (m): '))
        rho = float(input('Please specify the density of the fluid (kg/m^3): '))
        mu = float(input('Please specify the viscosity of the fluid (Pa*s): '))
        eps = float(input('Please specify the pipe surface roughness (m): '))
        L = float(input('Please specify the pipe length (m): '))
        deltaz = float(input('Please specify the elevation difference (m): '))
        Q = case2_SI(deltaP,D,rho,mu,eps,L,deltaz)
        print('Flow rate (m^3/s): ',Q)
    else :
        deltaP = float(input('Please specify the (-) delta pressure drop (psi): '))
        D = float(input('Please specify the pipe diameter (ft): '))
        rho = float(input('Please specify the density of the fluid (lbm/ft^3): '))
        mu = float(input('Please specify the viscosity of the fluid (lbm/(ft s)): '))
        eps = float(input('Please specify the pipe surface roughness (ft): '))
        L = float(input('Please specify the pipe length (ft): '))
        deltaz = float(input('Please specify the elevation difference (ft): '))
        Q = case2_FPS(deltaP,D,rho,mu,eps,L,deltaz)
        print('Flow rate (ft^3/s): ',Q)
else:
    if unit == 'SI':
        deltaP = float(input('Please specify the (-) delta pressure drop (Pa): '))
        Q = float(input('Please specify the flow rate(m^3/s): '))
        rho = float(input('Please specify the density of the fluid (kg/m^3): '))
        mu = float(input('Please specify the viscosity of the fluid (Pa*s): '))
        eps = float(input('Please specify the pipe surface roughness (m): '))
        L = float(input('Please specify the pipe length (m): '))
        deltaz = float(input('Please specify the elevation difference (m): '))
        D = case3_SI(deltaP,Q,rho,mu,eps,L,deltaz)
        print('Pipe diameter(m): ', D)
    else:
        deltaP = float(input('Please specify the (-) delta pressure drop (psi): '))
        Q = float(input('Please specify the flow rate(ft^3/s): '))
        rho = float(input('Please specify the density of the fluid (lbm/ft^3): '))
        mu = float(input('Please specify the viscosity of the fluid (lbm/(ft s)): '))
        eps = float(input('Please specify the pipe surface roughness (ft): '))
        L = float(input('Please specify the pipe length (ft): '))
        deltaz = float(input('Please specify the elevation difference (ft): '))
        D = case3_FPS(deltaP,Q,rho,mu,eps,L,deltaz)
        print('Pipe diameter(ft): ', D)