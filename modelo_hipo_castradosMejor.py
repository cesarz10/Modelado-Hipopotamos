##
from matplotlib import markers
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# para crecimiento poblacional sin control
def N_tn(I0, J0, A0, M0, control,n):  # nit = número de iteraciones/generaciones

    sacrificio = (
                control == 'S' )  # (booleano) checkea si el mecanismo de control introducido es sacrificio, mixto o emigración

    alpha = 0.5 * 0.8  # la proporcion es 1 hembra por cada 2 machos, y 4 adultos de cada 5 son reproductores
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

    I = np.eye(4, 4)  # matriz identidad
    Nt = np.vstack(np.array([I0, J0, A0, M0]))  # vector de distribución de edades
    N = np.sum(Nt)  # sumatoria de vector de distibución de edades

    Leslie = np.array([[pI, 0 , alpha * pA, 0], # matriz de Lefkovich #http://matema.ujaen.es/jnavas/web_recursos/archivos/matriciales/modelos%20matriciales%20tablas%20vida%20leslie.pdf
                       [aI, pJ,       0   , 0],
                       [ 0, aJ,      pA   , 0],
                       [ 0,  0,      aA   ,pM]])


    # si se usa el sacrificio, se saca n número de hipopótamos de todas las categorías por año (cte)
    if sacrificio:
        next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación
        arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número
        final_arr = []  # array al que se le va a append el # de individuos de cada categoría después de sacrificar

        for i in range(len(arr)):
            if i ==0:
                final_arr.append(arr[i])
            elif i == 3:
                final_arr.append(arr[i])

            else:
                final_arr.append(0) if arr[i] < n else final_arr.append(
                    arr[i] - n)  # se sacrifican n por generación, si hay menos de n entonces quedan 0 individuos

        return final_arr



        return final_arr

    next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación

    arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número
    return arr



# función para hallar valores de la siguiente generación teniendo en cuenta adultos castrados y no-castrados
def N_tn_cast(I0, J0, Ac0, A0, M0, control, b,n):  # nit = número de iteraciones/generaciones

    mixto = (
                control == 'M' )

    alpha = 1*0.5*0.8  # tasa de reproducción (1 por cada 0.5 mujeres y 0.8 adultos reprodutores)
    mI = 0.295 # prop NO supervivencia años 0-1
    mJA = 0.02 # prop NO supervivencia años 1-55

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
                       [ 0,       b*pJ  , sA,       0   , 0],
                       [ 0,       (1-b)*pJ,  0,      sA   , 0],
                       [ 0,        0    , pA,      pA   ,sM]])
    

    # si se usa el sacrificio, se saca n número de hipopótamos de todas las categorías por año (cte)
    if mixto:
        next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación
        arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número
        final_arr = []  # array al que se le va a append el # de individuos de cada categoría después de sacrificar
          # número de individuos a sacrificar por cada categoría (total sacrificados=n*2) | Adultos no castrados y juveniles

        for i in range(len(arr)):
            if i == 1 or i == 2 or i==3: # si es un adulto entonces se sacrifica
                final_arr.append(0) if arr[i] < n else final_arr.append(arr[i] - n)
            else: # si no es adulto no se sacrifica
                final_arr.append(arr[i])

        return final_arr


    next_it = Nt + ((K - N) / K) * np.dot((Leslie - I), (Nt))  # solo para la siguiente generación
    arr = np.around(next_it, 0)  # cogiendo el único vector de la matriz que da algún número

    return arr



# función para graficar el crecimiento poblacional
def graph(n_iter, I0, J0, A0, M0, control_met,b,n):
    if control_met =='C':
        Ac0 = 0
        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * (A0 + Ac0), np.ones(1) * M0
        Io, Jo, Aco, Ao, Mo = I0, J0, Ac0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None
        hipos = {}
        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = 'C'  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, Ac0, A0, M0 = N_tn_cast(I0, J0, Ac0, A0, M0, control,b,n)
            adultos = Ac0 + A0
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, adultos)
            M = np.append(M, M0)
            hipos["Inf"] = I
            hipos["Juv"] = J
            hipos["Cast"] = Ac0
            hipos["NoCast"] = A0
            hipos["Adu"] = adultos
            hipos["May"] = M


    elif control_met == 'M':
        Ac0 = 0
        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * (A0 + Ac0), np.ones(1) * M0
        Io, Jo, Aco, Ao, Mo = I0, J0, Ac0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None
        hipos = {}

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = 'M'  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, Ac0, A0, M0 = N_tn_cast(I0, J0, Ac0, A0, M0, control,b,n)
            adultos = Ac0 + A0
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, adultos)
            M = np.append(M, M0)
            hipos["Inf"] = I
            hipos["Juv"] = J
            hipos["Cast"] = Ac0
            hipos["NoCast"] = A0
            hipos["Adu"] = adultos
            hipos["May"] = M


    elif control_met == '':
        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * A0, np.ones(1) * M0
        Io, Jo, Ao, Mo = I0, J0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None
        hipos = {}

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = ''  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, A0, M0 = N_tn(I0, J0, A0, M0, control,n)
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, A0)
            M = np.append(M, M0)
            hipos["Inf"] = I
            hipos["Juv"] = J
            hipos["Adu"] = A0
            hipos["May"] = M

            
    elif control_met == 'S':

        I, J, A, M = np.ones(1) * I0, np.ones(1) * J0, np.ones(1) * A0, np.ones(1) * M0
        Io, Jo, Ao, Mo = I0, J0, A0, M0  # guadrando valores iniciales para título del gráfico
        control = None
        hipos = {}

        for i in range(n_iter):
            if i == 1:  # después de i generaciones se introduce el mecanismo de control
                control = 'S'  # !!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!! -> 'C' castración, 'S' sacrificio, 'M' mixto
            I0, J0, A0, M0 = N_tn(I0, J0, A0, M0, control,n)
            I = np.append(I, I0)
            J = np.append(J, J0)
            A = np.append(A, A0)
            M = np.append(M, M0)
            hipos["Inf"] = I
            hipos["Juv"] = J
            hipos["Adu"] = A0
            hipos["May"] = M
    print(hipos)
    return [I, J, A, M]



# valores iniciales de cada categoría
i, j, a, m = 2, 61, 39, 29

# resultados para los diferentes modelos de control
no_control = graph(80, i, j, a, m, '',0.25,8)
castracion = graph(50, i, j, a, m, 'C',0.25,8)
sacrificio = graph(50, i, j, a, m, 'S',0.25,8)
mixto = graph(50, i, j, a, m, 'M',0.4, 5) #Castran 0.25 Juveniles al año y se sacrifican 5 juveniles, y 5 adultos NO castrados


# subplot 2x2 de modelo sin control y combinación de modelos mixtos
figs, axs = plt.subplots(2, 2, figsize=(10,8))
figs.tight_layout(pad=3.5) # separando los subplots para que no queden tan encima
figs.subplots_adjust(top=0.9) # para que el título del subplot no se sobrelape con el de los axs

axs[0,0].plot(no_control[0], 'b', label="Infantes", marker='.', alpha=0.9)
axs[0,0].plot(no_control[1], 'r', label="Juveniles", marker='.', alpha=0.9)
axs[0,0].plot(no_control[2], 'g', label="Adultos", marker='.', alpha=0.9)
axs[0,0].plot(no_control[3], 'm', label="Mayores", marker='.', alpha=0.9)
axs[0,0].set_ylim(-10) # no me muestre nada por debajo del -10 en y
axs[0,0].set_xlim(-1) # no me muestre nada antes del -1 en x
axs[0,0].vlines(51, -10, 65, linewidth= 2,linestyle=':',color='black',alpha=0.4)
axs[0,0].hlines(65, 0, 51, linewidth= 2,linestyle=':',color='black',alpha=0.4)
axs[0,0].set_ylabel('Individuos')
axs[0,0].set_title('Sin control', y=1.012, fontsize=9)
axs[0,0].grid(alpha=0.3)


axs[0,1].plot(mixto[0], 'b', marker='^', markersize=4, alpha=1)
axs[0,1].plot(mixto[1], 'r', marker='^', markersize=4, alpha=1)
axs[0,1].plot(mixto[2], 'g', marker='^', markersize=4, alpha=1)
axs[0,1].plot(mixto[3], 'm', marker='^', markersize=4, alpha=1)
axs[0,1].plot(no_control[0], 'b', marker='.', alpha=0.2)
axs[0,1].plot(no_control[1], 'r', marker='.', alpha=0.2)
axs[0,1].plot(no_control[2], 'g', marker='.', alpha=0.2)
axs[0,1].plot(no_control[3], 'm', marker='.', alpha=0.2)
axs[0,1].set_xlim(0, 51)
axs[0,1].set_ylim(0, 90)
axs[0,1].set_title('Proporcion Castr=0.4 | Sacrificio/Emigración = 16 hipos', y=1.012, fontsize=9)
axs[0,1].grid(alpha=0.3)


mixto =  graph(50, i, j, a, m, 'M',0.3, 8) #BETA 0.5 y sacrificio 5
axs[1,0].plot(mixto[0], 'b', marker='^', markersize=4, alpha=1)
axs[1,0].plot(mixto[1], 'r', marker='^', markersize=4, alpha=1)
axs[1,0].plot(mixto[2], 'g', marker='^', markersize=4, alpha=1)
axs[1,0].plot(mixto[3], 'm', marker='^', markersize=4, alpha=1)
axs[1,0].plot(no_control[0], 'b', marker='.', alpha=0.2)
axs[1,0].plot(no_control[1], 'r', marker='.', alpha=0.2)
axs[1,0].plot(no_control[2], 'g', marker='.', alpha=0.2)
axs[1,0].plot(no_control[3], 'm', marker='.', alpha=0.2)
axs[1,0].set_xlim(0, 51)
axs[1,0].set_ylim(0, 90)
axs[1,0].set_xlabel('Años')
axs[1,0].set_ylabel('Individuos')
axs[1,0].set_title('Proporcion Castr=0.3 | Sacrificio/Emigración = 16 hipos', y=1.012, fontsize=9)
axs[1,0].grid(alpha=0.3)


mixto =  graph(50, i, j, a, m, 'M',0.2, 15) #BETA 0.5 y sacrificio 10
axs[1,1].plot(mixto[0], 'b', marker='^', markersize=4, alpha=1)
axs[1,1].plot(mixto[1], 'r', marker='^', markersize=4, alpha=1)
axs[1,1].plot(mixto[2], 'g', marker='^', markersize=4, alpha=1)
axs[1,1].plot(mixto[3], 'm', marker='^', markersize=4, alpha=1)
axs[1,1].plot(no_control[0], 'b', marker='.', alpha=0.2)
axs[1,1].plot(no_control[1], 'r', marker='.', alpha=0.2)
axs[1,1].plot(no_control[2], 'g', marker='.', alpha=0.2)
axs[1,1].plot(no_control[3], 'm', marker='.', alpha=0.2)
axs[1,1].set_xlim(0, 51)
axs[1,1].set_ylim(0, 90)
axs[1,1].set_xlabel('Años')
axs[1,1].set_title('Proporcion Castr=0.2 | Sacrificio/Emigración = 30 hipos', y=1.012, fontsize=9)
axs[1,1].grid(alpha=0.3)

# figs.suptitle(f'Crecimiento poblacional Hipopótamos  - I({int(i)}), J({int(j)}),  A({int(a)}), M({int(m)})', y=0.98)
figs.suptitle(f'Comparación de mecanismos de control mixtos', y=0.98)
figs.legend(loc='center right', bbox_to_anchor=(1.1, 0.6), fancybox=True, shadow=False, fontsize=9,framealpha=0.0)

plt.savefig(f'subplot_mixtos_Transp_I{int(i)}_J{int(j)}_A{int(a)}_M{int(m)}.png',
                transparent=True, bbox_inches='tight')
plt.show()


# # subplot 2x2 de modelo sin control, castracion, sacrificio, y mixto
# figs, axs = plt.subplots(2, 2, figsize=(10,8))
# figs.tight_layout(pad=3.5) # separando los subplots para que no queden tan encima
# figs.subplots_adjust(top=0.9) # para que el título del subplot no se sobrelape con el de los axs

# axs[0,0].plot(no_control[0], 'b', marker='.', alpha=0.9)
# axs[0,0].plot(no_control[1], 'r',marker='.', alpha=0.9)
# axs[0,0].plot(no_control[2], 'g',marker='.', alpha=0.9)
# axs[0,0].plot(no_control[3], 'm', marker='.', alpha=0.9)
# axs[0,0].set_ylim(-10) # no me muestre nada por debajo del -10 en y
# axs[0,0].set_xlim(-1) # no me muestre nada antes del -1 en x
# axs[0,0].vlines(51, -10, 65, linewidth= 2,linestyle=':',color='black',alpha=0.4)
# axs[0,0].hlines(65, 0, 51, linewidth= 2,linestyle=':',color='black',alpha=0.4)
# axs[0,0].set_ylabel('Individuos')
# axs[0,0].set_title('Sin control', y=1.012, fontsize=9)
# axs[0,0].grid(alpha=0.3)


# axs[0,1].plot(castracion[0], 'b', marker='s', markersize=3, alpha=0.9)
# axs[0,1].plot(castracion[1], 'r', marker='s', markersize=3, alpha=0.9)
# axs[0,1].plot(castracion[2], 'g', marker='s', markersize=3, alpha=0.9)
# axs[0,1].plot(castracion[3], 'm', marker='s', markersize=3, alpha=0.9)
# axs[0,1].plot(no_control[0], 'b', marker='.', alpha=0.2)
# axs[0,1].plot(no_control[1], 'r', marker='.', alpha=0.2)
# axs[0,1].plot(no_control[2], 'g', marker='.', alpha=0.2)
# axs[0,1].plot(no_control[3], 'm', marker='.', alpha=0.2)
# axs[0,1].set_xlim(0, 51)
# axs[0,1].set_ylim(0, 90)
# axs[0,1].set_title('Castración', y=1.012, fontsize=9)
# axs[0,1].grid(alpha=0.3)

# axs[1,0].plot(sacrificio[0], 'b', marker='^', markersize= 4,alpha=0.9)
# axs[1,0].plot(sacrificio[1], 'r', marker='^', markersize= 4,alpha=0.9)
# axs[1,0].plot(sacrificio[2], 'g', marker='^', markersize= 4,alpha=0.9)
# axs[1,0].plot(sacrificio[3], 'm', marker='^', markersize= 4,alpha=0.9)
# axs[1,0].plot(no_control[0], 'b', marker='.', alpha=0.2)
# axs[1,0].plot(no_control[1], 'r', marker='.', alpha=0.2)
# axs[1,0].plot(no_control[2], 'g', marker='.', alpha=0.2)
# axs[1,0].plot(no_control[3], 'm', marker='.', alpha=0.2)
# axs[1,0].set_xlim(0, 51)
# axs[1,0].set_ylim(0, 90)
# axs[1,0].set_xlabel('Años')
# axs[1,0].set_ylabel('Individuos')
# axs[1,0].set_title('Sacrificio', y=1.012, fontsize=9)
# axs[1,0].grid(alpha=0.3)

# axs[1,1].plot(castracion[0], 'b', marker='s', markersize=3, alpha=0.5)
# axs[1,1].plot(castracion[1], 'r', marker='s', markersize=3, alpha=0.5)
# axs[1,1].plot(castracion[2], 'g', marker='s', markersize=3, alpha=0.5)
# axs[1,1].plot(castracion[3], 'm', marker='s', markersize=3, alpha=0.5)
# axs[1,1].plot(sacrificio[0], 'b', marker='^',markersize=3,  alpha=0.5)
# axs[1,1].plot(sacrificio[1], 'r', marker='^', markersize=3, alpha=0.5)
# axs[1,1].plot(sacrificio[2], 'g', marker='^', markersize=3, alpha=0.5)
# axs[1,1].plot(sacrificio[3], 'm', marker='^', markersize=3, alpha=0.5)
# axs[1,1].plot(no_control[0], 'b', marker='.', alpha=0.1)
# axs[1,1].plot(no_control[1], 'r', marker='.', alpha=0.1)
# axs[1,1].plot(no_control[2], 'g', marker='.', alpha=0.1)
# axs[1,1].plot(no_control[3], 'm', marker='.', alpha=0.1)
# axs[1,1].set_xlim(0, 51)
# axs[1,1].set_ylim(0, 90)
# axs[1,1].set_xlabel('Años')
# axs[1,1].set_title('Comparación castración/sacrificio', y=1.012, fontsize=9)
# axs[1,1].grid(alpha=0.3)

# line1 = Line2D(range(1), range(1), color="black", marker='.', markersize=10) # sin control
# line2 = Line2D(range(1), range(1), color="black", marker='s') # castración
# line3 = Line2D(range(1), range(1), color="black", marker='^') # sacrificio
# line4 = Line2D(range(1), range(1), color="b") # infantes
# line5 = Line2D(range(1), range(1), color="r") # juveniles
# line6 = Line2D(range(1), range(1), color="g") # adultos
# line7 = Line2D(range(1), range(1), color="m") # mayores

# figs.suptitle(f'Comparación de los diferentes modelos de control')
# plt.legend((line4,line5,line6, line7, line1,line2,line3),('Infantes','Juveniles', 'Adultos', 'Mayores', 'Sin control','Castración', 'Sacrificio'), loc='center right', bbox_to_anchor=(1.4, 0.8), fancybox=True, shadow=False, fontsize=9,framealpha=0.0)
# # plt.legend((line4,line5,line6, line7),('Infantes','Juveniles', 'Adultos', 'Mayores'), loc='center right', bbox_to_anchor=(1.4, 0.3), fancybox=True, shadow=False, fontsize=9,framealpha=0.0)

# plt.savefig(f'subplot_modelos_Transp_I{int(i)}_J{int(j)}_A{int(a)}_M{int(m)}.png',
#                 transparent=True, bbox_inches='tight')
# plt.show()  