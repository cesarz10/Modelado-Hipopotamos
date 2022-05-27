import numpy as np
import matplotlib.pyplot as plt

# función para hallar valores de la siguiente generación
def N_tn(I0, J0, A0, M0, control, dict): # nit = número de iteraciones/generaciones

    castracion = (control == 'C' or control == 'M') # (booleano) checkea si el mecanismo de control introducido es castración
    sacrificio = (control == 'S' or control == 'M') # (booleano) checkea si el mecanismo de control introducido es sacrificio

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
        reproductores = (((A0*0.8) - dict['castrados']) / A0) if hipos_posibles <= dict['castrados'] else ((A0*0.8))
        print(reproductores)
        alpha = float(1*0.5*reproductores)


    I = np.eye(4, 4) # matriz identidad
    Nt = np.vstack(np.array([I0, J0, A0, M0])) # vector de distribución de edades
    N = np.sum(Nt) 

    Leslie = np.array([[  pI,    0, alpha*pA,  0], # matriz de Lefkovich #http://matema.ujaen.es/jnavas/web_recursos/archivos/matriciales/modelos%20matriciales%20tablas%20vida%20leslie.pdf
                       [  aI,    pJ,    0,     0],
                       [  0 ,    aJ,    pA,    0],
                       [  0 ,    0 ,    aA,   pM]])


    # si se usa el sacrificio, se saca n número de hipopótamos por año (cte)
    if sacrificio:
        next_it = Nt+((K-N) / K) * np.dot((Leslie- I), (Nt)) # solo para la siguiente generación
        arr = np.around(next_it,0) # cogiendo el único vector de la matriz que da algún número
        final_arr = [] # array al que se le va a append el # de individuos de cada categoría después de sacrificar
        n = 7 # número de individuos a sacrificar por cada categoría

        for i in range(len(arr)):
            final_arr.append(0) if arr[i] < n else final_arr.append(arr[i] - n) # se sacrifican n por generación, si hay menos de n entonces quedan 0 individuos
        
        return final_arr


    next_it = Nt+((K-N) / K) * np.dot((Leslie- I), (Nt))# solo para la siguiente generación
    arr = np.around(next_it,0)# cogiendo el único vector de la matriz que da algún número
    return arr


def graph(n_iter, I0, J0, A0, M0, dict2):
    I, J, A , M = np.ones(1)*I0, np.ones(1)*J0, np.ones(1)*A0, np.ones(1)*M0
    Io, Jo, Ao, Mo = I0, J0, A0, M0
    control = dict2['type']


    for i in range(n_iter):
        #if i == 40: # después de 40 generaciones se introduce el mecanismo de control
           # control = ''  # !!!!!!!!!!!! CAMBIAR A LA LETRA DEL MECANISMO DE CONTROL !!!!
        I0, J0, A0, M0 = N_tn(I0, J0, A0, M0, control, dict2['params'])
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
    plt.title(f'Crecimiento poblacional - {control} - I({int(Io)}), J({int(Jo)}), A({int(Ao)}), M({int(Mo)})')
    plt.xlabel('Generaciones')
    plt.ylabel('Individuos')
    hipos_posibles = dict2['params']['castrados']
    plt.savefig(f'Crecimiento_poblacional_control({control})_I{int(Io)}_J{int(Jo)}_A{int(Ao)}_M{int(Mo)}_{hipos_posibles}.png', transparent=True, bbox_inches='tight')
    plt.show()
presupuesto = 3622170000#CORNARE: https://www.cornare.gov.co/Acuerdos/Acuerdo_421_2021_cornare.pdf página 9
#presupuesto =492470000 #MIN AMBIENTE: https://www.minambiente.gov.co/wp-content/uploads/2022/04/Presupuesto-de-inversion-2022.pdf página 6
hipos_posibles = np.int32(presupuesto/25000000) #tambien puede ser 11000000, https://www.eltiempo.com/colombia/medellin/costo-de-exportar-los-hipopotamos-de-pablo-escobar-fuera-de-colombia-376448
print(hipos_posibles)

population = 133
# graph(100, population*0.02, population*0.46, population*0.30, population*0.22, {'type': 'C','params': {'castrados':hipos_posibles}}) #https://www.eltiempo.com/vida/medio-ambiente/hipopotamos-de-pablo-escobar-podrian-ser-declarados-especie-invasora-649043 #https://www.walshmedicalmedia.com/open-access/population-structure-of-the-common-hippopotamus-hippopotamus-amphibius-in-the-luangwa-river-zambia.pdf
graph(30, 2, 61, 39, 29, {'type': 'S', 'params': {'castrados':hipos_posibles}}) #https://www.eltiempo.com/vida/medio-ambiente/hipopotamos-de-pablo-escobar-podrian-ser-declarados-especie-invasora-649043 #https://www.walshmedicalmedia.com/open-access/population-structure-of-the-common-hippopotamus-hippopotamus-amphibius-in-the-luangwa-river-zambia.pdf
#graph(100, 0, 0, 4, 0)


# --------------------------------------------------------------------------------------------------------------