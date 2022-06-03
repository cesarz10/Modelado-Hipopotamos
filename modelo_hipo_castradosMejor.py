import numpy as np
import matplotlib.pyplot as plt

# para crecimiento poblacional sin control
def N_tn(I0, J0, A0, M0, control):  # nit = número de iteraciones/generaciones

    castracion = (
                control == 'C' or control == 'M')  # (booleano) checkea si el mecanismo de control introducido es castración o mixto
    sacrificio = (
                control == 'S' or control == 'M' or control == 'E')  # (booleano) checkea si el mecanismo de control introducido es sacrificio, mixto o emigración

    alpha = 1 * 0.5 * 0.8  # tasa de reproducción
    mI = 0.295  # tasa mort años 0-1
    mJ = 0.02  # tasa mort 1-55
    K = 1500  # capacidad de carga

    pI = 0  # porcentaje en grupo de edad
    pJ = 0.85
    pA = 0.95
    pM = 0.5

    aI = 1 - mI - pI  # porcentaje avance grupo de edad #http://www.conservationecologylab.com/uploads/1/9/7/6/19763887/lewison_2007.pdf
    aJ = 1 - mJ - pJ
    aA = 1 - mJ - pA

    # si se usa la castracion como mecanismo de control -> tasa de reproducción se ve reducida
    if castracion:
        alpha = alpha * 0.69  # castrando por proporción (69% capaz de reproducción)

    I = np.eye(4, 4)  # matriz identidad
    Nt = np.vstack(np.array([I0, J0, A0, M0]))  # vector de distribución de edades
    N = np.sum(Nt)  # sumatoria de vector de distibución de edades

    Leslie = np.array([[pI, 0 , alpha * pA, 0], # matriz de Lefkovich #http://matema.ujaen.es/jnavas/web_recursos/archivos/matriciales/modelos%20matriciales%20tablas%20vida%20leslie.pdf
                       [aI, pJ,       0   , 0],
                       [ 0, aJ,      pA   , 0],
                       [ 0,  0,      aA   ,pM]])

    # print(A0)

    # si se usa el sacrificio, se saca n número de hipopótamos de todas las categorías por año (cte)
    if sacrificio:
        next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación
        arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número
        final_arr = []  # array al que se le va a append el # de individuos de cada categoría después de sacrificar
        n = 4  # número de individuos a sacrificar por cada categoría (total sacrificados=n*4)

        for i in range(len(arr)):
            final_arr.append(0) if arr[i] < n else final_arr.append(
                arr[i] - n)  # se sacrifican n por generación, si hay menos de n entonces quedan 0 individuos

        return final_arr

    next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación

    arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número
    return arr



# función para hallar valores de la siguiente generación teniendo en cuenta adultos castrados y no-castrados
def N_tn_cast(I0, J0, Ac0, A0, M0, mixto):  # nit = número de iteraciones/generaciones

    sacrificio = mixto

    alpha = 1*0.5*0.8  # tasa de reproducción (1 por cada 0.5 mujeres y 0.8 adultos reprodutores)
    mI = 0.295 # prop NO supervivencia años 0-1
    mJA = 0.02 # prop NO supervivencia años 1-55
    b = 0.5 # proporción de juveniles castrados
    K = 1500  # capacidad de carga


    sI = 0   #proporcion que se queda en su grupo de edad
    sJ = 0.85
    sA = 0.95
    sM = 0.5

    pI = 1-mI-sI #porcentaje avance grupo de edad #http://www.conservationecologylab.com/uploads/1/9/7/6/19763887/lewison_2007.pdf
    pJ = 1-mJA-sJ
    pA = 1-mJA-sA
    pM = 1-mJA-sM
    # si se usa la castracion como mecanismo de control -> tasa de reproducción se ve reducida

    I = np.eye(5, 5)  # matriz identidad
    Nt = np.vstack(np.array([I0, J0, Ac0, A0, M0]))  # vector de distribución de edades
    N = np.sum(Nt)

    Leslie = np.array([[sI,        0    ,  0, alpha * pA, 0], # matriz de Lefkovich #http://matema.ujaen.es/jnavas/web_recursos/archivos/matriciales/modelos%20matriciales%20tablas%20vida%20leslie.pdf
                       [pI,       sJ    ,  0,       0   , 0],
                       [ 0,     b * pJ  , sA,       0   , 0],
                       [ 0, (1 - b) * pJ,  0,      sA   , 0],
                       [ 0,        0    , pA,      pA   ,sM]])
    

    # si se usa el sacrificio, se saca n número de hipopótamos de todas las categorías por año (cte)
    if sacrificio:
        next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación
        arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número
        final_arr = []  # array al que se le va a append el # de individuos de cada categoría después de sacrificar
        n = 4  # número de individuos a sacrificar por cada categoría (total sacrificados=n*4)

        for i in range(len(arr)):
            final_arr.append(0) if arr[i] < n else final_arr.append(arr[i] - n) 
            # se sacrifican n por generación, si hay menos de n entonces quedan 0 individuos

        return final_arr


    next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación
    arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número

    return arr



# función para graficar el crecimiento poblacional
def graph(n_iter, I0, J0, A0, M0, control_met):
    if control_met =='C':
        Ac0 = 0
        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * (A0 + Ac0), np.ones(1) * M0
        Io, Jo, Aco, Ao, Mo = I0, J0, Ac0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = 'C'  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, Ac0, A0, M0 = N_tn_cast(I0, J0, Ac0, A0, M0, mixto=False)
            adultos = Ac0 + A0
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, adultos)
            M = np.append(M, M0)

    elif control_met == 'M':
        Ac0 = 0
        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * (A0 + Ac0), np.ones(1) * M0
        Io, Jo, Aco, Ao, Mo = I0, J0, Ac0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = 'M'  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, Ac0, A0, M0 = N_tn_cast(I0, J0, Ac0, A0, M0, mixto=True)
            adultos = Ac0 + A0
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, adultos)
            M = np.append(M, M0)


    elif control_met == '':
        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * A0, np.ones(1) * M0
        Io, Jo, Ao, Mo = I0, J0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = control_met  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, A0, M0 = N_tn(I0, J0, A0, M0, control)
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, A0)
            M = np.append(M, M0)
            
    elif control_met == 'S':

        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * A0, np.ones(1) * M0
        Io, Jo, Ao, Mo = I0, J0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = control_met  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, A0, M0 = N_tn(I0, J0, A0, M0, control)
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, A0)
            M = np.append(M, M0)

    return [I, J, A, M]


# valores iniciales de cada categoría
i, j, a, m = 2, 61, 39, 29

# resultados para los diferentes modelos de control
no_control = graph(50, i, j, a, m, '')
castracion = graph(50, i, j, a, m, 'C')
sacrificio = graph(50, i, j, a, m, 'S')
mixto = graph(50, i, j, a, m, 'M')

# subplot 2x2 de modelo sin control, castracion, sacrificio, y mixto
figs, axs = plt.subplots(2, 2, figsize=(10,8), sharex=True)

axs[0,0].plot(no_control[0], 'b',label="Infantes", marker='^', alpha=0.8)
axs[0,0].plot(no_control[1], 'r',label="Juveniles", marker='^', alpha=0.8)
axs[0,0].plot(no_control[2], 'g',label="Adultos",marker='^', alpha=0.8)
axs[0,0].plot(no_control[3], 'm', label="Mayores",marker='^', alpha=0.8)
axs[0,0].set_xlabel('Años')
axs[0,0].set_ylabel('Individuos')
axs[0,0].set_title('Sin control')
axs[0,0].grid(alpha=0.3)


axs[0,1].plot(castracion[0], 'b', marker='.', alpha=0.9)
axs[0,1].plot(castracion[1], 'r', marker='.', alpha=0.9)
axs[0,1].plot(castracion[2], 'g', marker='.', alpha=0.9)
axs[0,1].plot(castracion[3], 'm', marker='.', alpha=0.9)
axs[0,1].plot(no_control[0], 'b', marker='^', alpha=0.1)
axs[0,1].plot(no_control[1], 'r', marker='^', alpha=0.1)
axs[0,1].plot(no_control[2], 'g', marker='^', alpha=0.1)
axs[0,1].plot(no_control[3], 'm', marker='^', alpha=0.1)
axs[0,1].set_xlabel('Años')
axs[0,1].set_ylabel('Individuos')
axs[0,1].set_title('Castración')
axs[0,1].grid(alpha=0.3)

axs[1,0].plot(sacrificio[0], 'b', marker='.', alpha=0.9)
axs[1,0].plot(sacrificio[1], 'r', marker='.', alpha=0.9)
axs[1,0].plot(sacrificio[2], 'g', marker='.', alpha=0.9)
axs[1,0].plot(sacrificio[3], 'm', marker='.', alpha=0.9)
axs[1,0].plot(no_control[0], 'b', marker='^', alpha=0.1)
axs[1,0].plot(no_control[1], 'r', marker='^', alpha=0.1)
axs[1,0].plot(no_control[2], 'g', marker='^', alpha=0.1)
axs[1,0].plot(no_control[3], 'm', marker='^', alpha=0.1)
axs[1,0].set_xlabel('Años')
axs[1,0].set_ylabel('Individuos')
axs[1,0].set_title('Sacrificio')
axs[1,0].grid(alpha=0.3)

axs[1,1].plot(mixto[0], 'b', marker='.', alpha=0.9)
axs[1,1].plot(mixto[1], 'r', marker='.', alpha=0.9)
axs[1,1].plot(mixto[2], 'g', marker='.', alpha=0.9)
axs[1,1].plot(mixto[3], 'm', marker='.', alpha=0.9)
axs[1,1].plot(no_control[0], 'b', marker='^', alpha=0.1)
axs[1,1].plot(no_control[1], 'r', marker='^', alpha=0.1)
axs[1,1].plot(no_control[2], 'g', marker='^', alpha=0.1)
axs[1,1].plot(no_control[3], 'm', marker='^', alpha=0.1)
axs[1,1].set_xlabel('Años')
axs[1,1].set_ylabel('Individuos')
axs[1,1].set_title('Mixto')
axs[1,1].grid(alpha=0.3)

figs.tight_layout()
figs.legend(loc='center', bbox_to_anchor=(1.1, 0.50), fancybox=True, shadow=False, fontsize=9,framealpha=0.0)

plt.savefig(f'subplot_completo_I{int(i)}_J{int(j)}_A{int(a)}_M{int(m)}.png',
                transparent=False, bbox_inches='tight')
plt.show()