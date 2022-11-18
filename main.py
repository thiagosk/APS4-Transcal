import numpy as np
from funcoesTermosol import importa
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')

x1 = N[0][0]
y1 = N[1][0]
x2 = N[0][1]
y2 = N[1][1]
x3 = N[0][2]
y3 = N[1][2]

L1 = ((x2 - x1)**2 + (y2-y1)**2)**0.5
L2 = ((x3 - x2)**2 + (y3-y2)**2)**0.5
L3 = ((x1 - x3)**2 + (y1-y3)**2)**0.5

sen1 = (y2-y1)/L1
cos1 = (x2-x1)/L1
sen2 = (y3-y2)/L2
cos2 = (x3-x2)/L2
sen3 = (y1-y3)/L3
cos3 = (x1-x3)/L3

K1 = Inc[0][2]*Inc[0][3]/L1*np.array([[cos1**2, cos1*sen1, -cos1**2, -cos1*sen1],
                      [cos1*sen1, sen1**2, -cos1*sen1, -sen1**2], 
                      [-cos1**2, -cos1*sen1, cos1**2, cos1*sen1], 
                      [-cos1*sen1, -sen1**2, cos1*sen1, sen1**2]])

K2 = Inc[1][2]*Inc[1][3]/L2*np.array([[cos2**2, cos2*sen2, -cos2**2, -cos2*sen2],
                      [cos2*sen2, sen2**2, -cos2*sen2, -sen2**2],
                      [-cos2**2, -cos2*sen2, cos2**2, cos2*sen2],
                      [-cos2*sen2, -sen2**2, cos2*sen2, sen2**2]])

K3 = Inc[2][2]*Inc[2][3]/L3*np.array([[cos3**2, cos3*sen3, -cos3**2, -cos3*sen3],
                      [cos3*sen3, sen3**2, -cos3*sen3, -sen3**2],
                      [-cos3**2, -cos3*sen3, cos3**2, cos3*sen3],
                      [-cos3*sen3, -sen3**2, cos3*sen3, sen3**2]])

def matriz_para_dicio(K, n):
    dicio = {}
    for linha in range(len(K)):
        for coluna in range(len(K[linha])):
            if n == 1:
                dicio[linha, coluna] = [(linha, coluna),K[linha][coluna]]
            if n == 2:
                dicio[linha, coluna] = [(linha+2, coluna+2),K[linha][coluna]]
            if n == 3:
                if linha == 0 or linha == 1:
                    if coluna == 0 or coluna == 1:
                        dicio[linha, coluna] = [(linha+4, coluna+4),K[linha][coluna]]
                    elif coluna == 2 or coluna == 3:
                        dicio[linha, coluna] = [(linha+4, coluna-2),K[linha][coluna]]
                elif linha == 2 or linha == 3:
                    if coluna == 0 or coluna == 1:
                        dicio[linha, coluna] = [(linha-2, coluna+4),K[linha][coluna]]
                    elif coluna == 2 or coluna == 3:
                        dicio[linha, coluna] = [(linha-2, coluna-2),K[linha][coluna]]
    return dicio

dicio1 = matriz_para_dicio(K1, 1)
dicio2 = matriz_para_dicio(K2, 2)
dicio3 = matriz_para_dicio(K3, 3)

Kg = np.zeros((6,6))
for linha in range(len(Kg)):
    for coluna in range(len(Kg[linha])):
        v = 0
        for valor in dicio1.values():
            if (linha, coluna) == valor[0]:
                v += valor[1]
                break
        for valor in dicio2.values():
            if (linha, coluna) == valor[0]:
                v += valor[1]
                break
        for valor in dicio3.values():
            if (linha, coluna) == valor[0]:
                v += valor[1]
                break
        Kg[linha][coluna] = v

Pg = []
for i in range(len(F)):
    if F[i][0] == 0:
        c = False
        if i <= nr:
            for j in R:
                if i == j[0]:
                    Pg.append("R")
                    c = True
        if c == False:
            Pg.append(0)
    else:
        Pg.append(F[i][0])