import math as m
import numpy as np

eps=10**(-4) #порядок точности
N=(25, 25) #кол-во разбиений
l=(1, 1) #стороны прямоугольной области
h=(l[0]/N[0],l[1]/N[1])

def k22(x, y):
    return 1 + x**2 + y**2

def ex_sol(x, y):
    return x*(x-1)*y*(y-1)

def f(x, y):
    return -2*y*(y-1)-2*x*(x-1)*(1-y+3*y**2+x**2)

def make_ex_sol ():
    u=np.zeros((N[0]+1, N[1]+1))
    for i in range (N[0]+1):
        for j in range(N[1]+1):
            u[i][j]=ex_sol(i*h[0], j*h[1])
    return u

def make_K22():
    K22=np.zeros((N[0]+1, N[1]+1))
    for i in range(N[0]+1):
        for j in range(N[1]+1):
            K22[i][j]=k22(i*h[0], j*h[1])
    return K22

def make_Au(u):
    Au_F=np.zeros((N[0]+1, N[1]+1))
    K22=make_K22()
    for i  in range(1, N[0]):
        for j in range(1, N[1]):
            L1=(u[i+1][j]-2*u[i][j]+u[i-1][j])/h[0]**2
            L2=(u[i][j+1]*(K22[i][j]+K22[i][j+1])/2-u[i][j]*(2*K22[i][j]+K22[i][j+1]+K22[i][j-1])/2+u[i][j-1]*(K22[i][j]+K22[i][j-1])/2)/h[1]**2
            Au_F[i][j] = -(L1 + L2)
    return Au_F

def make_F():
    F=np.zeros((N[0]+1, N[1]+1))
    for i  in range(1, N[0]):
        for j in range(1, N[1]):
            F[i][j]=f(i * h[0], j * h[1])
    return F
    
def make_B(F): 
    B=np.zeros((N[0]+1, N[1]+1))
    mu=[]
    for i in range(N[0]+1):
        mu.append([0]*(N[1]+1))
    for k1 in range(1,N[0]):
        for k2 in range(1, N[1]):
            for i in range(1, N[0]):    
                mu_k2=0
                for j in range(1, N[1]):
                    mu_k2+=F[i][j]*m.sin(k2*m.pi*j/N[1])
                mu[k1][k2]+=mu_k2*m.sin(k1*m.pi*i/N[0])
    for j in range(1, N[1]):
        for i in range(1, N[0]):
            for k2 in range(1, N[1]):
                v_k2=0
                for k1 in range(1, N[0]):
                    v_k2+=mu[k1][k2]/((4/h[0]**2)*(m.sin(k1*m.pi*h[0]/(2*l[0])))**2+4/(h[1]**2)*(m.sin(k2*m.pi*h[1]/(2.0 * l[1])))**2)*m.sin(k1*m.pi*i/N[0])
                B[i][j]+= 4*v_k2*m.sin(k2*m.pi*j/N[1])/(N[0]*N[1])
    return B

def norm(A):
    e=0
    for i in range(len(A)):
        for j in range(len(A[i])):
            if e<abs(A[i][j]):
                e=abs(A[i][j])
    return e

def make_alpha(r, y):
    Ay=make_Au(y)
    Ay_v=Ay.reshape(-1)
    r_v=r.reshape(-1)
    y_v=y.reshape(-1)
    alpha=np.dot(y_v,r_v)/np.dot(Ay_v,y_v)
    return alpha
    

def solve(): 
    e=10000
    alpha = 0.5
    err=[]
    u=np.zeros((N[0]+1, N[1]+1))
    count = 0
    while(e>eps):
        count+=1
        r=make_Au(u)-make_F()
        e = norm(r)
        err.append(e)
        y=make_B(r)
        alpha=make_alpha(r, y)
        u-=alpha*y

    for i in range(len(err)):
        print('||r||:')
        print(err[i])

    print('кол-во шагов: ', count)
    print('||v^n-v_ex||= ', norm(u-make_ex_sol())) 
        
    
solve()

