import math as m
import numpy as np
import matplotlib.pyplot as plt

a=0
b=1
#eps=10**(-2)

def f(x):#свободный член
    return x-1/2
    #return m.exp(-x)
    #return 5*x/6

def K(x, s): #ядро
    return 1
    #return x*m.exp(s)/2
    #return x*s/2

def make_U (h):
    R1=np.matrix([[h*K(i, a)/2] for i in np.arange(a, b+h,h)]) 
    R2=np.matrix([[2*h*K(i, j)/2 for j in np.arange(a+h, b,h)]
                 for i in np.arange(a, b+h,h)])
    R3=np.matrix([[h*K(i, b)/2] for i in np.arange(a, b+h,h)])
    R=np.hstack([np.hstack([R1, R2]), R3])
    n=int((b-a)/h)+1 #кол-во узлов
    E=np.eye(n)
    W=E-R
    print(W)
    F=np.matrix([[f(i)] for i in np.arange(a, b+h,h)])
    U=np.linalg.solve(W, F)
    return U

def u(x, h, U):
    u=h*K(x, a)*U[0]/2
    for i in range(int(a/h+1), int(b/h)):
        u+=2*h*K(x, i*h)*U[i]/2
    u+=h*K(x, b)*U[int(b/h)]/2
    u+=f(x)
    return float(u)

def delta_u(x, h):
    U1=make_U(h)
    U2=make_U(h/2)
    u1=u(x, h, U1)
    u2=u(x, h/2, U2)
    du=(float(u1-u2))**(2)
    return du

def gauss (func, a, b, h, n=5): 
    if n==3:
        t=(-(3/5)**(0.5), 0, (3/5)**(0.5))
        l=np.matrix([[2], [0], [2/3]])
    if n==4:
        sD=480**0.5
        t=(-(((30+sD)/70)**0.5), -(((30-sD)/70)**0.5), 
           (((30-sD)/70)**0.5), (((30+sD)/70)**0.5))
        l=np.matrix([[2], [0], [2/3], [0]])
    if n==5:
        sD=1120**0.5
        t=(-(((70+sD)/126)**0.5),-(((70-sD)/126)**0.5), 0,
           (((70-sD)/126)**0.5),(((70+sD)/126)**0.5))
        l=np.matrix([[2], [0], [2/3], [0], [2/5]])
    A=np.matrix([[t[i]**j for i in range(n)] for j in range(n)])
    c=np.linalg.inv(A)*l
    I=0
    for i in range(n):
        I+=c[i]*func((a+b)/2+(b-a)/2*t[i], h)
    return float(I*(b-a)/2)

def sleva_napravo(func, h1):
    I=0
    k=0
    eps_a=10**(-5)
    eps_o=10**(-5)
    p=2
    alpha=a
    betta=b
    while (alpha<betta):
       h=betta-alpha
       I_h=gauss(func, alpha, alpha+h, h1)
       I_h2=gauss(func, alpha, alpha+h/2, h1)+gauss(func, alpha+h/2, alpha+h, h1)
       delta_h=(I_h2-I_h)/((2**p)-1)
       k+=1
       if abs(delta_h) > max(eps_a, eps_o*abs(I_h2)):
           betta=(alpha+betta)/2
       else:
           I += I_h2+delta_h
           if abs(delta_h) > 1/(2**p)*max(eps_a, eps_o*abs(I_h2)):
               if betta+h > b:
                   alpha, betta = betta, b
               else:
                   alpha, betta = betta, betta+h
           else:
               if betta+2*h > b:
                   alpha, betta = betta, b
               else:
                   alpha, betta = betta, betta+2*h
    return I
def norma_c(h):
    U1=make_U(h)
    U2=make_U(h/2)
    lst=[abs(u(x, h, U1)-u(x, h/2, U2)) for x in np.arange(a, b, h)]
    delta_u=max(lst)
    return delta_u
def norma_l(h):
     return float(sleva_napravo(delta_u, h))**(0.5)
lst_h=[]
h=(b-a)/2
lst_h.append(h)
k=0
norm=norma_l(h)
while norm>10**(-5):
    print('h'+str(k)+' =',h)
    print('L2', (norm)**(0.5))
    print('C', norma_c(h))
    
    h=h/2
    lst_h.append(h)
    norm=norma_l(h)
    k+=1
print('h'+str(k)+' =',h)
print('L2',norm)
print('C', norma_c(h))
print('кол-во итераций =',k+1)
X=[x1 for x1 in np.arange(0, 1, 0.015625)]
y=[[u(x, t, make_U(t)) for x in X] for t in lst_h]
#y2=[x+m.exp(-x) for x in X]
y2=[x for x in X]
#y2=[x for x in X]
for t in range(len(lst_h)):
    text='u'+str(t)
    plt.plot(X, y[t],'--', label=text)
plt.plot(X, y2, label='fact')
plt.legend(fontsize=10)
plt.show()
