import numpy as np
from funcoesTermosol import importa
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada_exemplo.xlsx')


# Definindo elementos

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


# Calculando a matriz de rigidez de cada elemento

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


# Calculando a matriz de rigidez global 

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


# Formando o vetor global de forças concentradas

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


# Aplicando condição de contorno

lista_indices_nao_apoio_Pg = []
lista_indices_apoio_Pg = []
Pg_filtrado = []
for i in range(len(Pg)):
    if Pg[i] != "R":
        lista_indices_nao_apoio_Pg.append(i)
        Pg_filtrado.append(Pg[i])
    else:
        lista_indices_apoio_Pg.append(i)

Kg_filtrado_linha = []
for linha in range(len(Kg)):
    if linha in lista_indices_nao_apoio_Pg:
        Kg_filtrado_linha.append(Kg[linha])

Kg_filtrado = []
for linha in range(len(Kg_filtrado_linha)):
    l = []
    for coluna in range(len(Kg_filtrado_linha[linha])):
        if coluna in lista_indices_nao_apoio_Pg:
            l.append(Kg_filtrado_linha[linha][coluna])
    Kg_filtrado.append(l)


# Obtendo os deslocamentos nodais

multiplicador_das_incognitas = []
apoios = Pg_filtrado

contador = 0
while contador < len(Kg_filtrado):
    l = []
    for linha in range(len(Kg_filtrado)):
        for coluna in range(len(Kg_filtrado[linha])):
            if coluna == contador:
                l.append(Kg_filtrado[linha][coluna])
    multiplicador_das_incognitas.append(l)
    contador+=1

deslocamentos = np.linalg.solve(multiplicador_das_incognitas, apoios)
print(deslocamentos)


deslocamentos_expandido = np.zeros((len(Pg),1))
copia_deslocamentos = deslocamentos
for i in range(len(deslocamentos_expandido)):
    if i in lista_indices_nao_apoio_Pg:
        deslocamentos_expandido[i] = copia_deslocamentos[0]
        copia_deslocamentos = np.delete(copia_deslocamentos, 0)


# Obtendo as reações de apoio estrutural

apoios =  []
for linha in range(len(Kg)):
    calculo = [0]
    if linha in lista_indices_apoio_Pg:
        for coluna in range(len(Kg[linha])):
            calculo += Kg[linha][coluna] * deslocamentos_expandido[coluna]
        apoios.append([round(calculo[0],1)])


# Escrevendo o arquivo de saída

meuArquivo = open('saida.txt', 'w')
meuArquivo.write(f"Reacoes de apoio [N]\n{apoios}\n")
meuArquivo.write(f"\nDeslocamentos [m]\n{deslocamentos_expandido}\n")
meuArquivo.close()

