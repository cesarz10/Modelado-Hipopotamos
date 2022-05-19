import numpy as np
import matplotlib.pyplot as plt

def N_tn(I0, J0, A0, M0): # nit = número de iteraciones/generaciones

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
    aM = 1-mJ-pM

    I = np.eye(4, 4) # matriz identidad
    Nt = np.vstack(np.array([I0, J0, A0, M0])) # vector de distribución de edades
    N = np.sum(Nt) 

    Leslie = np.array([[  pI,    0, alpha*pA,  0], # matriz de Lefkovich #http://matema.ujaen.es/jnavas/web_recursos/archivos/matriciales/modelos%20matriciales%20tablas%20vida%20leslie.pdf
                       [  aI,    pJ,    0,     0],
                       [  0 ,    aJ,    pA,    0],
                       [  0 ,    0 ,    aA,   pM]])



    next_it = Nt+((K-N) / K) * np.dot((Leslie- I), (Nt))# solo para la siguiente generación

    # intentando para varias generaciones 
    # for i in range(nit):
    #     next_it = Nt + (((K-N) / K) * (Leslie - I) * (Nt))
    #     Nt = next_it[:,2].astype(np.int32)
    #     print(f'\ngeneracion {i}: \n{next_it}\n')

    arr = np.around(next_it,0)# cogiendo le único vector de la matriz que da algún número
    return arr


def graph(n_iter, I0, J0, A0, M0, arrP):
    I, J, A , M = np.ones(1)*2, np.ones(1)*13, np.ones(1)*82, np.ones(1)*3
    #figure = plt.figure()
    for i in range(n_iter):
        I0, J0, A0, M0 = N_tn(I0, J0, A0, M0)
        I = np.append(I, I0)
        J = np.append(J, J0)
        A = np.append(A, A0)
        M = np.append(M, M0)
    # plt.figure(figsize=((12,8)))
    plt.tight_layout()
    plt.plot(I, label = 'Infantes')
    plt.plot(J, label = 'Juveniles')
    plt.plot(A, label = 'Adultos')
    plt.plot(M, label = 'Ancianos')
    plt.legend(loc='center right', bbox_to_anchor=(1.30, 0.50), ncol=1, fancybox=True, shadow=False, fontsize=12, framealpha=0.2)
    plt.grid()
    plt.title('Crecimiento poblacional Hipopótamos')
    plt.xlabel('Generaciones')
    plt.ylabel('Individuos')
    plt.savefig(f'Crecimiento_poblacional_I{arrP[0]}, J{arrP[1]}, A{arrP[2]}, M{arrP[3]}.png', transparent=True, bbox_inches='tight')
    plt.show()

# graph(100, 2, 13, 82, 3)
arrP = [0, 0, 4, 0]
graph(50, 0, 0, 4, 0, arrP)

#Ejemplo Jensen
def N_tn(): # nit = número de iteraciones/generaciones
    K = 110  # capacidad de carga
    I = np.eye(4,4)

    #Nt = np.vstack(np.array([75, 57, 35, 22, 11, 4, 4]))  # vector de distribución de edades
    Nt = np.vstack(np.array([75, 57, 35, 22])) # vector de distribución de edades
    N = np.sum(Nt)
    '''
    Leslie = np.array([[  0.426    , 1.290, 1.296, 1.120, 1.126, 1.554, 0], # matriz de Leslie
                       [0.713, -1 , 0,0,0,0,0],
                       [ 0, 0.645, -1, 0,0,0,0],
                       [ 0,0,0.648,-1,0,0,0],
                       [0,0,0,0.560,-1,0,0],
                       [0,0,0,0,0.563,-1,0],
                       [0,0,0,0,0,0.777,-1]])

    '''
    Leslie = np.array([[0.426, 1.290, 1.296, 1.120], # matriz de Leslie
                       [0.713, -1 , 0,0],
                       [ 0, 0.645, -1, 0],
                       [ 0,0,0.648,-1]])

    next_it = Nt+((K-N) / K) * np.dot((Leslie-I), (Nt))# solo para la siguiente generación

    # intentando para varias generaciones
    # for i in range(nit):
    #     next_it = Nt + (((K-N) / K) * (Leslie - I) * (Nt))
    #     Nt = next_it[:,2].astype(np.int32)
    #     print(f'\ngeneracion {i}: \n{next_it}\n')

    arr = np.around(next_it,0)# cogiendo le único vector de la matriz que da algún número
    return arr

def graph(n_iter):
    #G1, G2, G3, G4, G5, G6, G7 = np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1)
    G1, G2, G3, G4 = np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1)
    '''
    for i in range(n_iter):
        G10, G20, G30, G40, G50, G60, G70 = N_tn()
        G1 = np.append(G1, G10)
        G2 = np.append(G2, G20)
        G3 = np.append(G3, G30)
        G4 = np.append(G4, G40)
        G5 = np.append(G5, G50)
        G6 = np.append(G6, G60)
        G7 = np.append(G7, G70)
    '''

    for i in range(n_iter):
        G10, G20, G30, G40 = N_tn()
        G1 = np.append(G1, G10)
        G2 = np.append(G2, G20)
        G3 = np.append(G3, G30)
        G4 = np.append(G4, G40)

    '''
    plt.plot(G1, label='1')
    plt.plot(G2, label='2')
    plt.plot(G3, label='3')
    plt.plot(G4, label='4')
    plt.plot(G5, label='5')
    plt.plot(G6, label='6')
    plt.plot(G7, label='7')
    '''
    plt.grid()
    plt.plot(G1, label='1')
    plt.plot(G2, label='2')
    plt.plot(G3, label='3')
    plt.plot(G4, label='4')
    plt.legend()
    plt.show()

#graph(25)