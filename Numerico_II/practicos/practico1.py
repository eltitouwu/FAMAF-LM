import numpy as np


#ej 1)

#a)
EPS=1e-9

def invalid(n):
    return np.array([np.nan]*n)

def sol_trinffil(A,b):
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float64)
    for i in range(m):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(np.inner(A[i,:i],x[:i])-b[i])>EPS):
                return invalid(m)
            continue
        x[i]=(b[i]-np.inner(A[i,:i+1],x[:i+1]))/A[i][i]
    return x

def sol_trinfcol(A,_b):
    b=np.ndarray.copy(_b)
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float64)
    for i in range(m):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(b[i])>EPS):
                return invalid(m)
            continue
        x[i]=b[i]/A[i][i]
        b[i+1:]-=np.inner(A[i+1:,i],x[i])
    return x

#b)

def sol_trsupfil(A,b):
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float64)
    for i in range(m-1,-1,-1):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(np.inner(A[i,i+1:],x[i+1:])-b[i])>EPS):
                return invalid(m)
            continue
        x[i]=(b[i]-np.inner(A[i,i+1:],x[i+1:]))/A[i][i]
    return x

def sol_trsupcol(A,_b):
    b=np.ndarray.copy(_b)
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float64)
    for i in range(m-1,-1,-1):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(b[i])>EPS):
                return invalid(m)
            continue
        x[i]=b[i]/A[i][i]
        b[:i]-=np.inner(A[:i,i],x[i])
    return x


#c)

from scipy.linalg import solve_triangular #para checkear

A=np.array([[1,0,0,0],[-1,1,0,0],[0,-1,1,0],[0,0,-2,2]],dtype=np.float64)
b=np.array([0,0,1,1],dtype=np.float64)
assert(np.array_equal(sol_trinfcol(A,b),solve_triangular(A,b,lower=True)))
assert(np.array_equal(sol_trinffil(A,b),solve_triangular(A,b,lower=True)))

A=np.array([[2,0,0,0],[-1,2,0,0],[3,1,-1,0],[4,1,-3,3]],dtype=np.float64)
b=np.array([2,3,2,9],dtype=np.float64)
assert(np.array_equal(sol_trinfcol(A,b),solve_triangular(A,b,lower=True)))
assert(np.array_equal(sol_trinffil(A,b),solve_triangular(A,b,lower=True)))

A=np.array([[9,2,4],[0,-6,3],[0,0,5]],dtype=np.float64)
b=np.array([18,-2,7],dtype=np.float64)
assert(np.array_equal(sol_trsupcol(A,b),solve_triangular(A,b,lower=False)))
assert(np.array_equal(sol_trsupfil(A,b),solve_triangular(A,b,lower=False)))

A=np.array([[1,2,-1,1],[0,1,0,-1],[0,0,-1,4],[0,0,0,1]],dtype=np.float64)
b=np.array([2,-1,0,0],dtype=np.float64)
assert(np.array_equal(sol_trsupcol(A,b),solve_triangular(A,b,lower=False)))
assert(np.array_equal(sol_trsupfil(A,b),solve_triangular(A,b,lower=False)))


#ej 5)

def cholesky(_X):
    X=np.ndarray.copy(_X)
    n=X.shape[0]
    assert(n==X.shape[1])
    G=np.zeros((n,n),dtype=X.dtype)
    for i in range(n):
        assert(X[i][i]>EPS)
        G[i][i]=np.sqrt(X[i][i])
        G[i,i+1:]=X[i,i+1:]/G[i][i]
        X[i+1:,i+1:]-=np.outer(G[i,i+1:],G[i,i+1:])
    return G

B=np.array([[4,-1,0],[-1,4,-1],[0,-1,4]],dtype=np.float64)
BB=np.block([[B,-np.eye(3),np.zeros((3,3))],[-np.eye(3),B,-np.eye(3)],[np.zeros((3,3)),-np.eye(3),B]])
G=cholesky(BB)


assert(np.allclose(G.T@G,BB,rtol=EPS,atol=EPS))


#ej 7)

def es_Sim(X):
    n=X.shape[0]
    if(n!=X.shape[1]): return False
    for i in range(n):
        for j in range(i):
            if(X[i,j]!=X[j,i]): return False
    return True
def es_SDP(_X):
    if(not es_Sim(_X)): return False
    X=np.ndarray.copy(_X)
    n=X.shape[0]
    
    G=np.zeros((n,n),dtype=X.dtype)
    for i in range(n):
        if(X[i,i]<=EPS): return False
        G[i,i]=np.sqrt(X[i,i])
        G[i,i+1:]=X[i,i+1:]/G[i,i]
        X[i+1:,i+1:]-=np.outer(G[i,i+1:],G[i,i+1:])
    return True

A=np.array([[4,0],[0,9]],dtype=np.float64)
assert(es_SDP(A))
G=cholesky(A)


#ej 8)

M=[]
L=['A','B','C','D']
A=np.array([[9,3,3],[3,10,5],[3,7,9]],dtype=np.float64)
M.append(A)
B=np.array([[4,2,6],[2,2,5],[6,5,29]],dtype=np.float64)
M.append(B)
C=np.array([[4,4,8],[4,-4,1],[8,1,6]],dtype=np.float64)
M.append(C)
D=np.array([[1,1,1],[1,2,2],[1,2,1]],dtype=np.float64)
M.append(D)
for i in range(4):
    if(es_SDP(M[i])):
        print(L[i]+" es SDP y su factor de Cholesky es:")
        print(cholesky(M[i]))
    else: print(L[i]+" no es SDP.")


#ej 9)

def sol_defpos(X,b):
    G=cholesky(X)
    return sol_trsupfil(G,sol_trinffil(G.T,b))



#ej 11)

L=np.float64(input("ingrese L para el ejercicio 11 (L>0): "))
assert(L>0.)
N=int(input("ingrese N para el ejercicio 11 (N>=2):"))
assert(N>=2)
A=np.zeros((N+2,N+2),dtype=np.float64)
b=np.array([0]*(N+2),dtype=np.float64)
h=L/(N+1)
A[0,0]=1
A[N+1,N+1]=1
b[0]=0
b[N+1]=0
for i in range(1,N+1):
    A[i,i-1]=-1/(h**2)
    A[i,i]=2/(h**2)
    A[i,i+1]=-1/(h**2)
    b[i]=(np.pi**2)*np.sin(np.pi*i*h)
A[1,0]=0
A[N,N+1]=0


v=(sol_defpos(A,b))

import matplotlib.pyplot as plt #por favor importar para graficar

x=[i*h for i in range(N+2)] #eje x.
y=[v[i] for i in range(N+2)] #eje y.
plt.plot(x,y) #grafica
plt.show()    #muestro la grafica




