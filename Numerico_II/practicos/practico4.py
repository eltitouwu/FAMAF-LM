import numpy as np
EPS=np.float128(10)**-15


#ej 3)
def givens_rot(a,b):
    if abs(b)<=EPS:
        return 1.,0.
    if abs(b)>abs(a):
        t=a/b
        s=-np.sign(b)/np.sqrt(1+t**2)
        c=s*t
    else:
        t=b/a
        c=np.sign(a)/np.sqrt(1+t**2)
        s=c*t
    return c,s

def qrgivens(A):
    M,N=A.shape
    Q=np.eye(M)
    R=np.ndarray.copy(A)
    P=min(N,M-1)
    for j in range(P):
        for k in range(j+1,M):
            if R[k,j]!=0:
                c,s=givens_rot(R[j,j],R[k,j])
                G=np.array([[c,s],[-s,c]])
                R[[j,k],j:]=G@R[[j,k],j:]
                R[k,j]=0.
                Q[:,[j,k]]=Q[:,[j,k]]@G.T
    return Q,R


def holder_ref(x):
    M=x.shape[0]
    s=0.
    for i in range(1,M): s+=x[i]**2
    if s<=EPS: return np.zeros((1,M)),0.
    m=np.sqrt(s+x[0]**2)
    if x[0]<=0: g=x[0]-m
    else: g=-s/(x[0]+m)
    r=2*g**2/(s+g**2)
    return np.block([np.eye(1),x[1:]/g]), r


def qrhholder(A):
    M,N=A.shape
    Q=np.eye(M)
    R=np.ndarray.copy(A)
    P=min(N,M)
    for j in range(P):
        u,r=holder_ref(R[j:,j])
        R[j:,j:]=R[j:,j:]-np.outer(r*u,u@R[j:,j:])
        Q[:,j:]=Q[:,j:]-np.outer(r*Q[:,j:]@u.T,u)
    return Q,R


A = np.array([[3.0,  5.0,  1.0],
              [4.0,  1.0,  2.0],
              [12.0, 2.0,  8.0]],dtype=np.float128)

Q,R=qrgivens(A)
# print(A)
# print(Q@R)
assert(np.allclose(Q@R,A,atol=EPS,rtol=EPS))
Q,R=qrhholder(A)
assert(np.allclose(Q@R,A,atol=EPS,rtol=EPS))


#ej 7)
from sol_tr import sol_trsupfil
def sol_cuadmin(A,b):
    Q,R=qrhholder(A.T@A)
    return sol_trsupfil(R,Q.T@(A.T@b))

#ej 8)

from sol_dlu import sol_dlu

from sol_tr import sol_defpos

for beta in range(1,10):
    eps=10**-beta
    A=np.array([[1.,1.],[eps,0.],[0.,eps]],dtype=np.float128)
    b=np.array([1.,1.,1.],dtype=np.float128)
    print('eps=',eps)
    AA=A.T@A
    Ab=A.T@b
    print(sol_dlu(AA,Ab),'LU')
    errLU=np.linalg.norm(A@sol_dlu(AA,Ab)-b,2)
    print(sol_defpos(AA,Ab),'Chol')
    errChol=np.linalg.norm(A@sol_defpos(AA,Ab)-b,2)
    print(sol_cuadmin(A,b),'QR')
    errQR=np.linalg.norm(A@sol_cuadmin(A,b)-b,2)
    errs=[(errLU,'LU'),(errChol,'Chol'),(errQR,'QR')]
    errs.sort()
    print('\n',errs,'\n')



#ej 9)
for n in [100,1000]:

    A=np.zeros((n,n-2))
    b=np.zeros(n)
    b[0]=1
    b[n-1]=1
    for i in range(n-2):
        A[i,i]=-1
        A[i+1,i]=2
        A[i+2,i]=-1
    AA=A.T@A
    Ab=A.T@b
    errChol=np.linalg.norm(A@sol_defpos(AA,Ab)-b,2)    
    errQR=np.linalg.norm(A@sol_cuadmin(A,b)-b,2)
    errs=[(errChol,'Chol'),(errQR,'QR')]
    errs.sort()
    print(errs)



#ej 11)
#a)
aux=np.pi*2./7.
Ac=np.array([[1,t,np.sin(aux*t),np.cos(aux*t)] for t in range(1,31)],dtype=np.float128)
b=np.loadtxt("y_mes1.txt",dtype=np.float128)

#b)
beta=sol_cuadmin(Ac,b)
print('\n',beta,'\n')

#c)
Ap=np.array([[1,t,np.sin(aux*t),np.cos(aux*t)] for t in range(31,61)],dtype=np.float128)
yp=Ap@beta

#d)
import matplotlib.pyplot as plt #por favor importar para graficar
y=np.loadtxt("y_mes2.txt",dtype=np.float128)

xAs=[t for t in range(1,31)] #eje x.
yAs=[b[t] for t in range(30)] #eje y.

plt.plot(xAs,yAs,color='black') #grafica


xBs=[t for t in range(31,61)] #eje x.
yBs=[y[t] for t in range(30)] #eje y.

plt.plot(xBs,yBs,color='blue',label='datos futuros') #grafica

xCs=[t for t in range(31,61)] #eje x.
yCs=[yp[t] for t in range(30)] #eje y.

plt.plot(xCs,yCs,color='red',label='datos predecidos') #grafica


plt.show()    #muestro la grafica

