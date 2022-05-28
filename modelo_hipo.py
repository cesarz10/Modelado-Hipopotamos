import numpy as np
import matplotlib.pyplot as plt

# función para hallar valores de la siguiente generación
def N_tn(I0, J0, A0, M0, control): # nit = número de iteraciones/generaciones

    castracion = (control == 'C' or control == 'M') # (booleano) checkea si el mecanismo de control introducido es castración o mixto
    sacrificio = (control == 'S' or control == 'M' or control == 'E') # (booleano) checkea si el mecanismo de control introducido es sacrificio, mixto o emigración

    alpha = 1*0.5*0.8  # tasa de reproducción
    mI = 0.295 # tasa mort años 0-1
    mJ = 0.02 # tasa mort 1-55
    K = 1500  # capacidad de carga

    pI = 0   #porcentaje en grupo de edad
    pJ = 0.85
    pA = 0.95
    pM = 0.5

    aI = 1-mI-pI #porcentaje avance grupo de edad #http://www.conservationecologylab.com/uploads/1/9/7/6/19763887/lewison_2007.pdf
    aJ = 1-mJ-pJ
    aA = 1-mJ-pA

    # si se usa la castracion como mecanismo de control -> tasa de reproducción se ve reducida
    if castracion:
        alpha = alpha * 0.69 # castrando por proporción (69% capaz de reproducción)

    I = np.eye(4, 4) # matriz identidad
    Nt = np.vstack(np.array([I0, J0, A0, M0])) # vector de distribución de edades
    N = np.sum(Nt) # sumatoria de vector de distibución de edades

    Leslie = np.array([[  pI,    0, alpha*pA,  0], # matriz de Lefkovich #http://matema.ujaen.es/jnavas/web_recursos/archivos/matriciales/modelos%20matriciales%20tablas%20vida%20leslie.pdf
                       [  aI,    pJ,    0,     0],
                       [  0 ,    aJ,    pA,    0],
                       [  0 ,    0 ,    aA,   pM]])

    # print(A0)

    # si se usa el sacrificio, se saca n número de hipopótamos de todas las categorías por año (cte)
    if sacrificio:
        next_it = Nt+((K-N) / K) * np.dot((Leslie- I), (Nt)) # solo para la siguiente generación
        arr = np.around(next_it,0) # cogiendo el único vector de la matriz que da algún número
        final_arr = [] # array al que se le va a append el # de individuos de cada categoría después de sacrificar
        n = 4 # número de individuos a sacrificar por cada categoría

        for i in range(len(arr)):
            final_arr.append(0) if arr[i] < n else final_arr.append(arr[i] - n) # se sacrifican n por generación, si hay menos de n entonces quedan 0 individuos
        
        return final_arr


    next_it = Nt+((K-N) / K) * np.dot((Leslie- I), (Nt)) # solo para la siguiente generación

    arr = np.around(next_it,0) # cogiendo el único vector de la matriz que da algún número
    return arr


# función para graficar el crecimiento poblacional
def graph(n_iter, I0, J0, A0, M0):
    I, J, A , M = np.ones(1)*I0, np.ones(1)*J0, np.ones(1)*A0, np.ones(1)*M0
    Io, Jo, Ao, Mo = I0, J0, A0, M0 # guadrando valores iniciales para título del gráfico
    control = None

    for i in range(n_iter):
        if i == 1: # después de i generaciones se introduce el mecanismo de control
            control = ''  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
        I0, J0, A0, M0 = N_tn(I0, J0, A0, M0, control)
        I = np.append(I, I0)
        J = np.append(J, J0)
        A = np.append(A, A0)
        M = np.append(M, M0)
        
    plt.tight_layout()
    plt.plot(I, label = 'Infantes', marker='.')
    plt.plot(J, label = 'Juveniles', marker='.')
    plt.plot(A, label = 'Adultos', marker='.')
    plt.plot(M, label = 'Ancianos', marker='.')
    plt.legend(loc='center right', bbox_to_anchor=(1.30, 0.50), ncol=1, fancybox=True, shadow=False, fontsize=12, framealpha=0.0)
    plt.grid()
    plt.title(f'Crecimiento poblacional - {control} - I({int(Io)}), J({int(Jo)}), A({int(Ao)}), M({int(Mo)})') # título con mecanismo de control
    # plt.title(f'Crecimiento poblacional - I({int(Io)}), J({int(Jo)}), A({int(Ao)}), M({int(Mo)})') # título de crecimiento poblacional sin mecanismo de control
    plt.xlabel('Años')
    plt.ylabel('Individuos')
    plt.savefig(f'Crecimiento_poblacional_control({control})_I{int(Io)}_J{int(Jo)}_A{int(Ao)}_M{int(Mo)}.png', transparent=True, bbox_inches='tight')
    plt.show()



graph(90, 2, 61, 39, 29) # -> 39 * 0.31 = 12.09 
# graph(90, 1371, 61, 39, 29) # Equilibrio
# graph(30, 0, 0, 4, 0)



# --------------------------------------------------------------------------------------------------------------