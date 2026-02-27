from math import * #por favor importar, en caso contrario usar serie_seno del primer parcial, y redefinir sin usando mi serie_seno con n=10, o definirse un poly con los coeficientes a gusto de la función sin. También redefinir la constante pi.
import copy #innecesario para el parcial de hoy, pero me quedo más tranquilo sabiendo que los poly copian sus coeficientes en vez de agarrarlos por referencia de otro poly u otra lista.

class poly:
    
    def __init__(self, coef): #coef es la lista de coeficientes de menor a mayor grado!!!!
        self.c=copy.deepcopy(coef)
    
    def __add__(p1, p2):
        r=poly([0. for i in range(max(p1.len(),p2.len()))])
        for i in range(r.len()):
            if(p1.len()>i): r.c[i]+=p1.c[i]
            if(p2.len()>i): r.c[i]+=p2.c[i]
        while r.len()>1 and abs(r.c[-1])==0.0: 
            r.c.pop() #Elimino coeficientes de mayor grado que se cancelaron entre p1 y p2, salvo que el polinomio sea el polinomio 0.
        return r
    
    def __mul__(p1,p2):
        r=poly([0. for i in range(p1.len()+p2.len()-1)])
        for i in range(p1.len()):
            for j in range(p2.len()): r.c[i+j]+=p1.c[i]*p2.c[j] #impresionante que la suma se codea más largo que el producto jaja.
        return r
    def len(p1):
        return len(p1.c) 
    def horn(p1,x):
        r=0.
        for i in range(p1.len()-1,-1,-1): r=r*x+p1.c[i] 
        return r
    
    def to_str(p1):
        s=""
        for i in range(p1.len()):
            s+="("+str(p1.c[i])+")"
            if(i>=1): s+=" x"
            if(i>1): s+="^"+str(i) 
            if(i<p1.len()-1): s+=" +"
        return s
    
    def deriv(p1):
        if(p1.len()==1): return poly([0.])
        return poly([p1.c[i]*i for i in range(1,p1.len())])
    def Integ(p1):
        P1=poly([0.])
        P1.c+=[p1.c[i]/(i+1) for i in range(p1.len())]
        return P1
#Ejercicio 1 =============================================================================

def simpson(fun, a, b, n):
    r=0.
    h=(b-a)/(2*n)
    x=a
    for i in range(n):
        r+=fun(x)+4*fun(x+h)+fun(x+2*h)
        x+=2*h
    return r*h/3

def n_pol(p1,a,b): #explicado en instrucciones2.txt.
    if(p1.len()<4): return 1
    p2=poly(p1.c)
    for i in range(p2.len()): p2.c[i]=abs(p2.c[i])
    for i in range(4): p2=p2.deriv()
    N=(1/2) * ((((b-a)**5)*1e5/180)*p2.horn(max(abs(a),abs(b))))**(1/4) #error de simpson.
    return int(max(1,ceil(N))) #no quiero N=0.

def n_sen(p1,a,b): #explicado en instrucciones2.txt.
    p2=poly([0]); auxp=poly(p1.c)
    for i in range(auxp.len()): auxp.c[i]=abs(auxp.c[i])
    Comb4=[1,4,6,4,1]
    for i in range(5): 
        p2=p2+(auxp*poly([Comb4[i]]))
        auxp=auxp.deriv()
    N=(1/2) * ((((b-a)**5)*1e5/180)*p2.horn(max(abs(a),abs(b))))**(1/4) #error de simpson.
    return int(max(1,ceil(N))) #no quiero N=0.

def POL_BOG(w,a,b,n):
    def prodvec(p1,p2):
        auxf=lambda x: (w.horn(x))*(p1.horn(x))*(p2.horn(x))
        return simpson(auxf,a,b,n_pol(w*p1*p2,a,b))
    r=[poly([1.])]                                      #phi_0
    if(n==0): return r      #terminé.
    x=poly([0,1])
    Bi=prodvec(x*r[0],r[0])/prodvec(r[0],r[0])
    r.append(poly([-Bi,1.]))                            #phi_1
    if(n==1): return r      #terminé.
    for k in range(2,n+1):
        Bi=prodvec(x*r[-1],r[-1])/prodvec(r[-1],r[-1])  
        Ci=prodvec(x*r[-1],r[-2])/prodvec(r[-2],r[-2])
        phik=(poly([-Bi,1.])*r[-1])+(poly([-Ci])*r[-2]) #phi_k   #esta operación es O(n) pues el producto de polinomios es de un polinomio de grado a lo sumo n con un polinomio de grado 1, cuya complejidad es O(n).
        r.append(phik)
    return r

BOG=POL_BOG(poly([1.]), 0, pi/2, 2)
coef=[simpson(lambda x:sin(x)*phi.horn(x),0,pi/2,n_sen(phi,0,pi/2))/simpson(lambda x: phi.horn(x)**2,0,pi/2,n_pol(phi*phi,0,pi/2)) for phi in BOG] #resuelvo el sistema diagonal de una, pues este ejercicio no me piden que aplique ningún método especial.

cuad_min_ej1=poly([0.])
for i in range(len(BOG)):
    cuad_min_ej1=cuad_min_ej1+(poly([coef[i]])*BOG[i]) #reconstruyo la rta.

print("EJ1: el polinomio cuadratico que aproxima a f(x)=sen(x) en el intervalo [0,pi/2] en el sentido de cuadrados minimos es:")
print("p(x)="+cuad_min_ej1.to_str(),"\n")

#Ejercicio 2 =============================================================================

def Lin_Tr(A,b,inf):
    n=len(A)
    U=n if inf else -1
    L=0 if inf else n-1
    D=1 if inf else -1 
    x=[0. for i in range(n)]
    for i in range(L,U,D):
        aux=0.
        for j in range(L,i,D):
            aux+=A[i][j]*x[j]
        x[i]=(b[i]-aux)/A[i][i]
    return x

def LU_fact(A):
    n=len(A)
    U=[[0. for j in range(n)] for i in range(n)]
    L=[[0. for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i,n):
            S=0.
            for k in range(0,i):
                S+=L[i][k]*U[k][j]
            U[i][j]=A[i][j]-S
        L[i][i]=1
        for j in range(i+1,n): 
            S=0.
            for k in range(0,i):
                S+=L[j][k]*U[k][i]
            L[j][i]=(A[j][i]-S)/U[i][i]
    return (L,U)

def Lin_LU(A,B):
    n=len(A)
    q=len(B) #q de queries
    X=[]
    (L,U)=LU_fact(A)
    for i in range(q):
        y=Lin_Tr(L,B[i],True)
        X.append(Lin_Tr(U,y,False))
    return X
#BOG=[poly([float(i==j) for j in range(i+1)]) for i in range(3)] #descomentar si quiere que A no sea diagonal.
b=[simpson(lambda x:sin(x)*phi.horn(x),0,pi/2,n_sen(phi,0,pi/2)) for phi in BOG] 
A=[[simpson(lambda x:BOG[i].horn(x)*BOG[j].horn(x),0,pi/2,n_pol(BOG[i]*BOG[j],0,pi/2)) for j in range(len(BOG))] for i in range(len(BOG))] #Ahora sí armo la matriz de las ecuaciones        
x=Lin_LU(A,[b])[0] #resuelvo el sistema lineal con LU.
cuad_min_ej2=poly([0.])
for i in range(len(BOG)):
    cuad_min_ej2=cuad_min_ej2+(poly([x[i]])*BOG[i]) #reconstruyo la respuesta.

print("EJ2: el polinomio cuadratico que aproxima a f(x)=sen(x) en el intervalo [0,pi/2] en el sentido de cuadrados minimos es:")
print("p(x)="+cuad_min_ej2.to_str(),"\n")

#Ejercicio 3 =============================================================================

def Ab_splincub(_a,_b,n,fun,nat,d2fun=None):
    x=[_a+i*(_b-_a)/(n) for i in range(n+1)]
    A=[[0. for j in range(4*n)] for i in range(4*n)]
    b=[0. for j in range(4*n)]
    for i in range(n): 
        aux1=1.
        aux2=1.
        for j in range(4): #Condiciones de interpolación y C0
            A[4*i][j+4*i]=aux1  
            A[4*i+1][j+4*i]=aux2 
            aux1*=x[i]
            aux2*=x[i+1]
        b[4*i]=fun(x[i])
        b[4*i+1]=fun(x[i+1])
        if(i==n-1): break    
        aux1=1.
        for j in range(1,4): #Condiciones de C1 y C2
            A[4*i+2][4*i+j]=aux1*j               #C1
            A[4*i+2][4*(i+1)+j]=-aux1*j          #C1
            A[4*i+3][4*i+j]=aux1*j*(j-1)         #C2
            A[4*i+3][4*(i+1)+j]=-aux1*j*(j-1)    #C2
            aux1*=x[i+1]
    A[4*n-2][2]=2.; A[4*n-2][3]=3*2*_a;          
    A[4*n-1][4*n-2]=2.; A[4*n-1][4*n-1]=3*2*_b; 
    if not nat: b[4*n-2]=d2fun(_a); b[4*n-1]=d2fun(_b)
    return (A,b)

class splincub:
    def __init__(self,_a,_b,n,fun,nat,d2fun=None):  #O((40n)^3)!!!   #n es cantidad de polinomios!!
        self.nat=nat
        self._a=_a
        self._b=_b
        self.n=n
        self.xi=[_a+_a+i*(_b-_a)/(n) for i in range(n+1)]
        (A,b)=Ab_splincub(_a,_b,n,fun,nat,d2fun)
        x=Lin_LU(A,[b])[0]
        self.pol=[]
        for i in range(n):
            ls=[]
            for j in range(4):
                ls.append(x[4*i+j])
            self.pol.append(poly(ls))
    def eval(splin,x):
        if(x<splin._a or x>splin._b): assert(0)
        izq=0; der=splin.n-1
        while(izq<=der):
            medio=(izq+der)//2
            if(x>=splin.xi[medio]): izq=medio+1
            else: der=medio-1
        return splin.pol[der].horn(x)
    def to_str(splin):
        r=""
        for i in range(splin.n):
            r+=splin.pol[i].to_str()
            r+=" para x en ["+str(splin.xi[i])+","+str(splin.xi[i+1])+"]\n"
        r+="Spline "+("natural" if splin.nat else "correcto")+"\n"
        return r

Splinej3=splincub(0,pi/2,4,sin,False,lambda x: -sin(x))
print("EJ3: el Spline Cubico que aproxima a f(x)=sen(x) en el intervalo [0,pi/2] con 5 nodos es:")
print("S(x)=\n"+Splinej3.to_str()+"\n")

#Ejercicio 4 ========================================================================================
import matplotlib.pyplot as plt #Por favor importar para poder graficar!!!!   #esta línea es la que traba el código jaja.

fondo_negro="Y"
fondo_negro=str(input("¿Quiere el gráfico con fondo negro (fachero y recomendado)? Y/n:"))
if(fondo_negro!="n"):
    face, ax = plt.subplots(facecolor="black")
    ax.set_facecolor('black')
    for ch in ('x','y'): plt.tick_params(axis=ch, colors='white')
    for bord in ('top', 'bottom', 'left', 'right'): ax.spines[bord].set_color('white')

n_pts=int(1e3)
x1=[i*pi/(2*(n_pts-1)) for i in range(n_pts)]
y1=[cuad_min_ej2.horn(v) for v in x1]
y2=[Splinej3.eval(v) for v in x1]
y3=[sin(v) for v in x1]

plt.plot(x1,y1,c="limegreen",label="Cuadrados Minimos")
plt.plot(x1,y2,c="blue",label="Spline")
plt.plot(x1,y3,c="red",linestyle=(5, (10, 5)),label="Seno")
plt.grid(True)
plt.legend()
plt.show()

#Extra ================================================================
Fact=[1] #lista de factoriales: posicion -> posicion!.
def init_fact(maxv):
    for i in range(len(Fact),maxv): #completo mi lista de factoriales.
        Fact.append(Fact[-1]*i)  #Fact[-1] -> ultimo elemento de Fact.
    return 

def Integ_pol(p1,a,b):
    P1=p1.Integ()
    return P1.horn(b)-P1.horn(a)

def Integ_senpol(p1,a,b): #La antiderivada es de la pinta P1S*sen + P1C*cos
    P1S=poly([0. for i in range(p1.len())]) 
    P1C=poly(P1S.c)
    if(p1.len()>len(Fact)): init_fact(p1.len())
    for i in range(p1.len()):
        for j in range(i+1):
            if j&1: P1S.c[i-j]+=p1.c[i]*(Fact[i]/Fact[i-j])*(-1 if (j>>1)&1 else 1) # + o - sen dependiendo de la congruencia de j mod 4
            else:   P1C.c[i-j]+=p1.c[i]*(Fact[i]/Fact[i-j])*(1 if (j>>1)&1 else -1) # - o + cos dependiendo de la congruencia de j mod 4
    return (P1S.horn(b)*sin(b)-P1S.horn(a)*sin(a))+(P1C.horn(b)*cos(b)-P1C.horn(a)*cos(a))


def POL_BOG2(w,a,b,n):
    def prodvec(p1,p2): return Integ_pol(w*p1*p2,a,b)
    r=[poly([1.])]                                      #phi_0
    if(n==0): return r      #terminé.
    x=poly([0,1])
    Bi=prodvec(x*r[0],r[0])/prodvec(r[0],r[0])
    r.append(poly([-Bi,1.]))                            #phi_1
    if(n==1): return r      #terminé.
    for k in range(2,n+1):
        Bi=prodvec(x*r[-1],r[-1])/prodvec(r[-1],r[-1])  
        Ci=prodvec(x*r[-1],r[-2])/prodvec(r[-2],r[-2])
        phik=(poly([-Bi,1.])*r[-1])+(poly([-Ci])*r[-2]) #phi_k   #esta operación es O(n) pues el producto de polinomios es de un polinomio de grado a lo sumo n con un polinomio de grado 1, cuya complejidad es O(n).
        r.append(phik)
    return r
Quiero_Extra="N"
Quiero_Extra=str(input("¿Quiere el extra? y/N:"))
if(Quiero_Extra=="y"):
    BOG=POL_BOG2(poly([1.]), 0, pi/2, 2)
    coef=[Integ_senpol(phi,0,pi/2)/Integ_pol(phi*phi,0,pi/2) for phi in BOG] #resuelvo el sistema diagonal de una, pues este ejercicio no me piden que aplique ningún método especial.
    cuad_min_ext=poly([0.])
    for i in range(len(BOG)):
        cuad_min_ext=cuad_min_ext+(poly([coef[i]])*BOG[i]) #reconstruyo la rta.

    print("EXTRA: el polinomio cuadratico que aproxima a f(x)=sen(x) en el intervalo [0,pi/2] en el sentido de cuadrados minimos es:")
    print("p(x)="+cuad_min_ext.to_str(),"\n")
    y4=[cuad_min_ext.horn(v) for v in x1]
    y1=[cuad_min_ej1.horn(v) for v in x1]
    if(fondo_negro!="n"):
        face, ax = plt.subplots(facecolor="black")
        ax.set_facecolor('black')
        for ch in ('x','y'): plt.tick_params(axis=ch, colors='white')
        for bord in ('top', 'bottom', 'left', 'right'): ax.spines[bord].set_color('white')

    plt.plot(x1,y1,c="limegreen",label="Cuadrados Minimos")
    plt.plot(x1,y2,c="blue",label="Spline")
    plt.plot(x1,y3,c="red",linestyle=(5, (10, 5)),label="Seno")
    plt.plot(x1,y4,c="purple",linestyle=(5, (10, 5)), label="Cuadrados Minimos Extra")
    plt.grid(True)
    plt.legend()
    plt.show()