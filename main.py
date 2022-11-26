import numpy as np
from funcoesTermosol import importa
from numpy.linalg import det
from matplotlib import pyplot as plt
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')

no_indice = {}

for i in range(1, nn + 1):
    no_indice[i] = (2*i-2,2*i-1)

class elemento:
    def __init__(self, posicao, no0, nof, E, A):
        self.posicao = posicao
        self.no0 = no_indice[no0]
        self.nof = no_indice[nof]
        self.x0 = N[0][int(no0-1)]
        self.y0 = N[1][int(no0-1)]
        self.xf = N[0][int(nof-1)]
        self.yf = N[1][int(nof-1)]
        self.L = ((self.xf - self.x0)**2 + (self.yf-self.y0)**2)**0.5
        self.sen = (self.yf-self.y0)/self.L
        self.cos = (self.xf-self.x0)/self.L
        self.E = E
        self.A = A
        self.K = E*A/self.L*np.array([[self.cos**2, self.cos*self.sen, -self.cos**2, -self.cos*self.sen],
                                      [self.cos*self.sen, self.sen**2, -self.cos*self.sen, -self.sen**2], 
                                      [-self.cos**2, -self.cos*self.sen, self.cos**2, self.cos*self.sen], 
                                      [-self.cos*self.sen, -self.sen**2, self.cos*self.sen, self.sen**2]])
        self.tensao = 0
        self.deformacao = 0
        self.forca = 0


def criacao_elementos():
    elementos = []  
    for i in range(nm):
        e = elemento(i+1, Inc[i][0], Inc[i][1], Inc[i][2], Inc[i][3])
        elementos.append(e)
    return elementos

elementos = criacao_elementos()

def matriz_para_dicio(elementos):
    for element in elementos:
        a = element.no0[0], element.no0[1], element.nof[0], element.nof[1]
        k_indice = {}
        for i, j in zip(element.K, a):
            for c, b in zip(i, a):
                k_indice[(j, b)] = c
        element.dicio = k_indice

matriz_para_dicio(elementos)

def calcula_Kg(elementos):
    tam = no_indice[nn][1] + 1
    Kg = np.zeros((tam,tam))
    
    for element in elementos:
        for i in range(tam):
            for j in range(tam):
                indice = (i,j)
                if indice in element.dicio:
                    Kg[i][j] += element.dicio[indice]
    return Kg

kg = calcula_Kg(elementos)

def calcula_Pg_total():
    Pg = []
    for i in range(len(F)):
        if F[i][0] == 0:
            c = False
            for j in R:
                if i == j[0]:
                    Pg.append("R")
                    c = True
            if c == False:
                Pg.append(0)
        else:
            Pg.append(F[i][0])
    return Pg

Pg_total = calcula_Pg_total()

def kg_filtrado(kg, Pg_total):
    indices = [i for i, j in enumerate(Pg_total) if j == "R"]
    tam = (no_indice[nn][1] + 1 - len(indices))
    Kg_novo = np.zeros((tam,tam))

    a = 0
    b = 0
    for i, linha in enumerate(kg):
        if i not in indices:
            for j, coluna in enumerate(linha):
                if j not in indices:
                    Kg_novo[a][b] = coluna
                    b += 1
            a += 1
            b = 0
    
    Pg = [i for i in Pg_total if i != "R"]
                
    return Kg_novo, Pg
    
    
kg_fil, Pg_fil = kg_filtrado(kg, Pg_total)

def gauss_seidel(Pg_fil, kg_fil, tol, max_iter):
    indices = [i for i, j in enumerate(Pg_total) if j == "R"]
    x = np.zeros(kg_fil.shape[0])
    
    erros = []
    
    for i in range(max_iter):
        for j in range(kg_fil.shape[0]):
            x[j] = (Pg_fil[j] - np.dot(kg_fil[j, :j], x[:j]) - np.dot(kg_fil[j, j+1:], x[j+1:]))/kg_fil[j, j]
        erro = np.linalg.norm(np.dot(kg_fil, x) - Pg_fil)
        erro = erro/np.linalg.norm(Pg_fil)
        erros.append(erro)
        if erro < tol:
            break
    
    x_expandido = x.copy()
    for i in indices:
        x_expandido = np.insert(x_expandido, i, 0)
    
    return x, x_expandido, erros
    
x, x_expandido, erros = gauss_seidel(Pg_fil, kg_fil, 1e-6, 1000)

def calcula_reacoes(x_expandido, kg):
    indices = [i for i, j in enumerate(Pg_total) if j == "R"]
    reacoes = np.dot(kg, x_expandido)
    reacoes = [reacoes[i] for i in indices]
    return reacoes

reacoes = calcula_reacoes(x_expandido, kg)

def tensão_deformação(elementos, x_expandido):
    for element in elementos:
        _cos = element.cos
        _sen = element.sen
        a = element.no0[0], element.no0[1], element.nof[0], element.nof[1]
        u = np.array([x_expandido[i] for i in a])
        element.tensao = (element.E/element.L)*np.dot(np.array([-_cos, -_sen, _cos, _sen]), u)
        element.deformacao = (1/element.L)*np.dot(np.array([-_cos, -_sen, _cos, _sen]), u)
        element.forca = element.tensao*element.A
    
tensão_deformação(elementos, x_expandido)

tensao = [element.tensao for element in elementos]
forca = [element.forca for element in elementos]
deformacao = [element.deformacao for element in elementos]

print("Reacoes de apoio [N]")
print(reacoes)

print("\n")
print("Deslocamentos [m]")
print(x_expandido)

print("\n")
print("Deformacoes []")
print(deformacao)

print("\n")
print("Forcas internas [N]")
print(forca)

print("\n")
print("Tensoes internas [Pa]")
print(tensao)

testex = N[0] + x_expandido[::2]
testey = N[1] + x_expandido[1::2]

def plota(N,Inc):
    # Numero de membros
    nm = len(Inc[:,0])
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt

#    plt.show()
    fig = plt.figure()
    # Passa por todos os membros
    for i in range(nm):
        
        # encontra no inicial [n1] e final [n2] 
        n1 = int(Inc[i,0])
        n2 = int(Inc[i,1])        

        plt.plot([N[0,n1-1],N[0,n2-1]],[N[1,n1-1],N[1,n2-1]],color='r',linewidth=3)


    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.grid(True)
    plt.axis('equal')

deforma = [testex, testey]
deforma = np.array(deforma)

plota(N, Inc)

plota(deforma, Inc)

plt.show()
