
eps=10**(-6)
a=0
b=1
y0=1
count = 0

def f(x, y):
    global count
    count+=1
    return y
    
def find_y (h, x, y):
    k1 = h*f(x, y)
    k2 = h*f(x + h/2, y + k1/2)
    k3 = h*f(x + h, y - k1 + 2*k2)
    return y+(k1+4*k2+k3)/6

def sleva_narpavo ():
    x=a
    y=y0
    h=b-a
    while (x!=b):
        y_nxt=find_y(h, x, y)
        y_sm=find_y(h/2, x, y)
        y_nxt_sm=find_y(h/2, x+h/2, y_sm)
        p=2
        rho=abs((y_nxt-y_nxt_sm)/(1-2**(-p)))
        if rho<eps*h/(b-a):
            y=y_nxt
            if x+h<b:
                x+=h
                h*=2
            else:
                h=b-x
                x=b 
             
        else:
            h=h/2 
    return y

print('y(b)=', sleva_narpavo())
print('кол-во вызовов:',count)
