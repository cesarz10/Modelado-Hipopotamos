##
import numpy as np
import matplotlib.pyplot as plt


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
def N_tn_cast(I0, J0, Ac0, A0, M0):  # nit = número de iteraciones/generaciones


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
    

    next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación
    arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número

    return arr


# función para graficar el crecimiento poblacional
def graph(n_iter, I0, J0, A0, M0,control_met):
    if control_met =='C':
        Ac0 = 0
        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * (A0 + Ac0), np.ones(1) * M0
        Io, Jo, Aco, Ao, Mo = I0, J0, Ac0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = 'C'  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, Ac0, A0, M0 = N_tn_cast(I0, J0, Ac0, A0, M0)
            adultos = Ac0 + A0
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, adultos)
            M = np.append(M, M0)
            marcador = '.'

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
            marcador = 's'
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
            marcador = '^'

    plt.tight_layout()

    plt.plot(I, 'b',label="Infantes", marker=marcador, alpha=0.35)
    plt.plot(J, 'r',label="Juveniles", marker=marcador, alpha=0.35)
    plt.plot(A, 'g',label="Adultos",marker=marcador, alpha=0.35)
    plt.plot(M, 'm', label="Mayores",marker=marcador, alpha=0.35)
    plt.legend(loc='center right', bbox_to_anchor=(1.30, 0.50), ncol=1, fancybox=True, shadow=False, fontsize=9,
               framealpha=0.0)
    plt.grid()
    plt.title(
        f'Crecimiento poblacional - {control} - I({int(Io)}), J({int(Jo)}),  A({int(Ao)}), M({int(Mo)})')  # título con mecanismo de control
    # plt.title(f'Crecimiento poblacional - I({int(Io)}), J({int(Jo)}), A({int(Ao)}), M({int(Mo)})') # título de crecimiento poblacional sin mecanismo de control
    plt.xlabel('Años')
    plt.ylabel('Individuos')
    plt.savefig(f'Crecimiento_poblacional_control({control})_I{int(Io)}_J{int(Jo)}_A{int(Ao)}_M{int(Mo)}.png',
                transparent=False, bbox_inches='tight')
    plt.show()


graph(80, 2, 61, 39, 29, '')  # -> 39 * 0.31 = 12.0





"""
def mostrar_menu() -> None:
    print("\nMenú de opciones:")
    print("1. Elegir parametros iniciales")
    print("2. Elegir metodo de control")
    print("3. Elija el numero de años para la simulación")
    print("4. Ver grafica")
    print("5. Salir")


def iniciar_aplicacion() -> None:
    print("Inicializando aplicacion")
    continuar = True
    metodo = ''
    while continuar:
        mostrar_menu()
        opcion = int(input("Ingrese la opcion que desea ejecutar: "))
        if opcion == 1:
            Infantes = int(input("Ingrese el numero de infantes iniciales: "))
            Juveniles = int(input("Ingrese el numero de Juveniles iniciales: "))
            Adultos = int(input("Ingrese el numero de Adultos iniciales: "))
            Mayores = int(input("Ingrese el numero de Mayores iniciales: "))
        elif opcion == 2:
            metodo = input("Ingrese el metodo de control adecuado: Castracion->(C) | Sin control->() | Sacrificio->(S)").strip()
        elif opcion == 3:
            años = int(input("Ingrese un numero de años para simular: "))
        elif opcion == 4:
            print(metodo)

            graph(años, Infantes, Juveniles, Adultos, Mayores, metodo)
        elif opcion == 5:
            continuar = False
        else:
            print("Ingresó una opción no válida, por favor elija una de las opciones del menú")


iniciar_aplicacion()

"""