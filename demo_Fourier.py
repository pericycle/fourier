import numpy as np
import matplotlib.pyplot as plt
import time


# animate the Fourier approximation
def plotAnimFourier(case,L,n,delay):
    # plots the function and all its Fourier approx
    # such that 0<=j<=n
    for j in np.arange(n+1):
        plotFourier(case,L,j)
        plt.title(str(j))
        plt.show()
        time.sleep(delay)

# plot one Fourier approximation of order n
def plotFourier(case,L,n):
    # plots the function f and n-term Fourier approx
    a, b = fourierCoeffs(case, L, n)
    npts = 1000
    x = np.linspace(-L,L,npts)
    f = funcgen(case)
    y = np.zeros((x.shape))
    
    for j in np.arange(x.size):
        y[j]=f(x[j])  
    ymin = y.min()
    ymax = y.max()
    xmin = x.min()
    xmax = x.max()
    offsetx = 0.2*(xmax-xmin)
    offsety = 0.2*(ymax-ymin)
    y1 = ymin-offsety
    y2 = ymax+offsety
    x1 = xmin-offsetx
    x2 = xmax+offsetx
    fsum = a[0]*np.ones((x.shape))
    for j in np.arange(1,a.size):
        fsum += a[j]*np.cos(j*x*np.pi/L)+ b[j]*np.sin(j*x*np.pi/L)
    plt.plot(x[:500],y[:500],linewidth='3',label='Function',color='blue')
    plt.plot(x[500:],y[500:],linewidth='3',label='', color='blue')
    plt.ylim(y1,y2)
    plt.xlim(x1,x2)
    plt.plot(x,fsum,color='red',linewidth='1.5',label='Fourier series')
    plt.legend()

# compute the Fourier coefficients using Simpson's quadrature
def fourierCoeffs(case,L,n):
    f = funcgen(case)
    eps = 1E-14
    npts = 1001
    x1 = np.linspace(-L,-eps,npts)
    x2 = np.linspace(eps,L,npts)
    x = np.hstack([x1,x2]).flatten()
    y = np.zeros((x.shape)).flatten()
    for j in np.arange(2*npts):
        y[j] = f(x[j])
    u = np.pi*x/L
    cons = np.ones((n+1,1))/L
    cons[0]=0.5*cons[0]
    A = np.zeros((n+1,1)).flatten()
    B = np.copy(A)
    for j in np.arange(n+1):
        fsin = y*np.sin(j*u)
        fcos = y*np.cos(j*u)
        B[j] = cons[j]*(simpson(fsin[:npts],-L,-eps,npts) + \
                     simpson(fsin[npts:],eps,L,npts))
        A[j] = cons[j]*(simpson(fcos[:npts],-L,-eps,npts) + \
                     simpson(fcos[npts:],eps,L,npts))
    return A, B
        

def simpson(f,a,b,n):
    # Simpson's rule to compute Fourier coefficients
    # n is odd
    if not isinstance(f,np.ndarray):
        y = np.zeros((n,1)).flatten()
        x = np.linspace(a,b,n)
        for j in np.arange(n):
            y[j] = f(x[j])
    else:
        y = f
        
    h = (b-a)/(n-1);
    integral = h*(y[0]+y[n-1])/3;
    k = 1;
    while (k < n-1):
        if (k%2==0):
            integral += 2*h*y[k]/3;
        else:
            integral += 4*h*y[k]/3;
        k+=1;
    return integral

def funcgen(case):
    if case=='jump':
        def f(x):
            if x>0:
                return 1
            else:
                return -1
        return f
    elif case=='line':
        def f(x):
            return x
        return f
    elif case=='parabola':
        def f(x):
            return x**2
        return f
    elif case=='even':
        def f(x):
            return np.exp(-np.abs(x))
        return f
    elif case=='sin_abs':
        def f(x):
            return np.sin(10*x)*np.abs(x)
        return f
    elif case=='sin_decay':
        def f(x):
            if x>0:
                return 5*np.sin(10*x)*np.exp(-1.2*x)
            else:
                return 0
        return f
    elif case=='sin_exp':
        def f(x):
            if x<0:
                return -2*np.sin(5*x)*np.exp(0.15*x)
            else:
                return x**2*np.sin(10*x)-8
        return f
    elif case=='sin':
        def f(x):
            return np.sin(x)
        return f
    elif case=='cos13':
        def f(x):
            return np.cos(13*x)
        return f
    elif case=='wiggly':
        def f(x):
            func = 2*np.exp(np.abs(x))*np.sin(5*x)**2
            if x<0:
                return func
            else:
                return -func-20
        return f
    elif case=='sin_ext':
        def f(x):
            if x<0:
                return 0
            else:
                return np.sin(x)
        return f
    elif case=='parajump':
        def f(x):
            if x<0:
                return (x+2)**2
            else:
                return 8-x
        return f
    elif case=='abs':
        def f(x):
            return np.abs(x)
        return f
    


        
        
            
            
            
    
    
    
