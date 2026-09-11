import numpy as np
import matplotlib.pyplot as plt #por favor importar para graficar
#ej 10) b)
N=10000
ang=[i*2*np.pi/N for i in range(N)]

for eps in [0.25,0.125,0.0625,1e-5]:
    vs=[np.array([np.cos(th),np.sin(th)]) for th in ang]
    xs=[v[0] for v in vs] #eje x.
    ys=[v[1] for v in vs] #eje y.
    ax=plt.subplot()
    ax.set_aspect('equal')
    plt.plot(xs,ys,color='black') #grafica
    A=np.array([[1,1-eps],[0,1]])

    xAs=[(A@v)[0] for v in vs] #eje x.
    yAs=[(A@v)[1] for v in vs] #eje y.

    plt.plot(xAs,yAs,color='red') #grafica


    B=np.array([[1/eps,0],[0,eps]])

    xBs=[(B@v)[0] for v in vs] #eje x.
    yBs=[(B@v)[1] for v in vs] #eje y.

    plt.plot(xBs,yBs,'blue') #grafica
    
    plt.show()    #muestro la grafica


#ej) 11)

A=np.loadtxt("A_dataset.txt",dtype=np.float64)
b=np.loadtxt("b_dataset.txt",dtype=np.float64)

from sol_egauss import sol_egauss
x=sol_egauss(A,b)
E=np.random.rand(A.shape[0],A.shape[1])

xs=[]
ys=[]
for beta in range(1,11):
    eps=1/beta
    B=A+eps*E
    y=sol_egauss(B,b)
    dx=np.linalg.norm(y-x,2)/np.linalg.norm(x,2)
    dA=np.linalg.norm(B-A,2)/np.linalg.norm(A,2)
    xs.append(dA)
    ys.append(dx)

plt.plot(xs,ys)
plt.xscale("log")
plt.yscale("log")
plt.grid(True, which="both", ls="--")
plt.show()
