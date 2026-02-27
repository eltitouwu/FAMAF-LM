

#ejercicio 1. ==========================================================================================================================================================

Fact=[1] #lista de factoriales: posicion -> posicion!.
def init_fact(maxv):
    for i in range(len(Fact),maxv): #completo mi lista de factoriales.
        Fact.append(Fact[-1]*i)  #Fact[-1] -> ultimo elemento de Fact.
    return 


def serie_seno(x):
    r=0. #inicializo el resultado de la funcion.
    n=5 #primeros n terminos del polinomio de taylor centrado en 0 de sen(x) (fijado n=5).

    if len(Fact)<2*n: #calculo los factoriales que me falten.
        init_fact(2*n)

    x2=x*x
    for i in range(n-1,-1,-1):  #algoritmo de Horn, modificado por particularidades de este polinomio.
        r*=x2
        if i&1 : # es lo mismo que i (mod 2), y recordar que 0 es lo mismo que False, y todo int distinto de 0 se considera True
            r-=1/Fact[2*i+1]   #(-1^impar)
        else :
            r+=1/Fact[2*i+1]   #(-1^par)
    r*=x
    return r 


#ejercicio 2. ==========================================================================================================================================================

import matplotlib.pyplot as plt #por favor importar para graficar

x=[i/100 for i in range(641)] #eje x.
y=[serie_seno(v) for v in x] #eje y.
plt.plot(x,y) #grafica
plt.show()    #muestro la grafica

#ejercicio 3. ==========================================================================================================================================================

def rsecante(fun,x0,x1,err,mit):

    hx=[x1]
    f0=fun(x0)
    f1=fun(x1)
    hf=[f1]

    for it in range(mit):

        if abs(f1)<err :
            break 

        x0=x0-(f0*(x0-x1))/(f0-f1) #guardo el nuevo x_i en el valor viejo de x_{i-2}.
        f0=fun(x0) #recalculo su evaluacion en la funcion.

        x0,x1 = x1,x0
        f0,f1 = f1,f0 #intercambio posiciones.

        hx.append(x1)
        hf.append(f1)

    return (hx,hf)


#ejercicio 4. ==========================================================================================================================================================

def rnewton(fun,x0,err,mit): 
    #me defino el metodo de Newton dados: 
    #   (1)fun es una funcion fun:x -> (f(x),f'(x)) 
    #   (2)x0 es el punto inicial del metodo. 
    #   (3)err es el error que busco en el rango de iteraciones. 
    #   (4)mit es la maxima iteracion del metodo.
    #debido a que lo necesito para busqueda_ceros.
    
    hx=[x0] 
    f0=fun(x0)
    hf=[f0[0]]
    #las tres lineas anteriores son un poco innecesarias pero las dejo para que los metodos sean similares.

    for it in range(mit):

        if abs(f0[0])<err:
            break

        x0=x0-f0[0]/f0[1] #calculo el nuevo punto.
        f0=fun(x0) #recalculo su evaluacion en la funcion.

        hx.append(x0)
        hf.append(f0[0])

    return (hx,hf)

def busqueda_ceros(fun,x0,x1,err,mit): #asumo que fun me devuelve (f(x),f'(x)).
    
    #Secante.
    def newfun(x):
        return fun(x)[0]
    (hx_sec,hf_sec)=rsecante(newfun,x0,x1,err,mit)
    print("Metodo de la secante encontro a la raiz ",hx_sec[-1]," en ",len(hx_sec)-1," iteraciones")
    
    #Newton.
    (hx_newton,hf_newton)=rnewton(fun,x0,err,mit)
    print("Metodo de Newton encontro a la raiz ",hx_newton[-1]," en ",len(hx_newton)-1," iteraciones")
    
    if abs(hf_sec[-1])<abs(hf_newton[-1]): #devuelvo la mejor raiz.
        return hx_sec[-1]
    else :
        return hx_newton[-1]



#ejercicio 5. ==========================================================================================================================================================

def serie_coseno(x): #me lo defino para tener la derivada de serie_seno
    r=0. #inicializo el resultado de la funcion.
    n=5 #primeros n terminos del polinomio de taylor centrado en 0 de cos(x) (fijado n=5).

    if len(Fact)<2*n: #calculo los factoriales que me falten.
        init_fact(2*n)

    x=x*x
    for i in range(n-1,-1,-1):  # algoritmo de Horn, modificado por particularidades de este polinomio.
        r*=x
        if i&1 : # es lo mismo que i (mod 2) .
            r-=1/Fact[2*i]
        else :
            r+=1/Fact[2*i]
    return r 

def sen_cos(x): #me lo defino para el Newton, pues sen'=cos
    return (serie_seno(x),serie_coseno(x))

busqueda_ceros(sen_cos,3,6,1e-5,100)
#requiere 5 pasos para el metodo de la secante.
#requiere 2 pasos para el metodo de newton.

busqueda_ceros(sen_cos,4.5,6,1e-5,100)
#requiere 7 pasos para el metodo de la secante.
#requiere 6 pasos para el metodo de newton.
#ninguno de los dos devuelven una aproximacion para una raiz de sen(x).
#esto es debido a que las series de taylor de sen y cos aproximan mal en ese rango 
#con los primeros 5 terminos no nulos.


#extra. ==========================================================================================================================================================

#desactivar los comentarios de las funciones busqueda_ceros2 llamadas abajo del todo
def deriv(fun,x): #aproximo la derivada de f en x
	return (fun(x*(1+1e-8))-fun(x))/(x*1e-8)


def busqueda_ceros2(fun,x0,x1,err,mit): #asumo que fun me devuelve f(x).
    
    #Secante.
    (hx_sec,hf_sec)=rsecante(fun,x0,x1,err,mit)
    print("Metodo de la secante encontro a la raiz ",hx_sec[-1]," en ",len(hx_sec)-1," iteraciones")
    
    #Newton.
    def newfun(x):
    	return (fun(x),deriv(fun,x))
    (hx_newton,hf_newton)=rnewton(newfun,x0,err,mit)
    print("Metodo de Newton encontro a la raiz ",hx_newton[-1]," en ",len(hx_newton)-1," iteraciones")
    
    if abs(hf_sec[-1])<abs(hf_newton[-1]): #devuelvo la mejor raiz.
        return hx_sec[-1]
    else :
        return hx_newton[-1]

#busqueda_ceros2(serie_seno,3,6,1e-5,100)


#busqueda_ceros2(serie_seno,4.5,6,1e-5,100)


