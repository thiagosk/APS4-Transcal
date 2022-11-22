import numpy as np
from funcoesTermosol import importa
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')


# Definindo elementos

def matriz_para_dicio(K, posicao):
    dicio = {}
    for linha in range(len(K)):
        for coluna in range(len(K[linha])):
            if posicao == 1:
                dicio[linha, coluna] = [(linha, coluna),K[linha][coluna]]
            if posicao == 2:
                dicio[linha, coluna] = [(linha+2, coluna+2),K[linha][coluna]]
            if posicao == 3:
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

class elemento:
    def __init__(self, posicao, no0, nof, E, A):
        self.posicao = posicao
        self.no0 = no0
        self.nof = nof
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
        self.dicio = matriz_para_dicio(self.K, posicao)

elementos = []
for i in range(nm):
    e = elemento(i+1, Inc[i][0], Inc[i][1], Inc[i][2], Inc[i][3])
    elementos.append(e)

trelicas = []
for element in elementos:
    trelica = []
    for e2 in elementos:
        if e2.no0 == element.nof and e2.no0 != element.no0:
            element2 = e2
            break
    for e3 in elementos:
        if e3.no0 == element2.nof and e3.no0 != e2.no0:
            element3 = e3
            break
    trelica.append(element.posicao)
    trelica.append(element2.posicao)
    trelica.append(element3.posicao)
    if sorted(trelica) not in trelicas:
        trelicas.append(trelica)

print(trelicas)


# x1 = N[0][0]
# y1 = N[1][0]
# x2 = N[0][1]
# y2 = N[1][1]
# x3 = N[0][2]
# y3 = N[1][2]

# L1 = ((x2 - x1)**2 + (y2-y1)**2)**0.5
# L2 = ((x3 - x2)**2 + (y3-y2)**2)**0.5
# L3 = ((x1 - x3)**2 + (y1-y3)**2)**0.5

# sen1 = (y2-y1)/L1
# cos1 = (x2-x1)/L1
# sen2 = (y3-y2)/L2
# cos2 = (x3-x2)/L2
# sen3 = (y1-y3)/L3
# cos3 = (x1-x3)/L3


# Calculando a matriz de rigidez de cada elemento

# K1 = Inc[0][2]*Inc[0][3]/L1*np.array([[cos1**2, cos1*sen1, -cos1**2, -cos1*sen1],
#                                       [cos1*sen1, sen1**2, -cos1*sen1, -sen1**2], 
#                                       [-cos1**2, -cos1*sen1, cos1**2, cos1*sen1], 
#                                       [-cos1*sen1, -sen1**2, cos1*sen1, sen1**2]])

# K2 = Inc[1][2]*Inc[1][3]/L2*np.array([[cos2**2, cos2*sen2, -cos2**2, -cos2*sen2],
#                                       [cos2*sen2, sen2**2, -cos2*sen2, -sen2**2],
#                                       [-cos2**2, -cos2*sen2, cos2**2, cos2*sen2],
#                                       [-cos2*sen2, -sen2**2, cos2*sen2, sen2**2]])

# K3 = Inc[2][2]*Inc[2][3]/L3*np.array([[cos3**2, cos3*sen3, -cos3**2, -cos3*sen3],
#                                       [cos3*sen3, sen3**2, -cos3*sen3, -sen3**2],
#                                       [-cos3**2, -cos3*sen3, cos3**2, cos3*sen3],
#                                       [-cos3*sen3, -sen3**2, cos3*sen3, sen3**2]])


# Calculando a matriz de rigidez global 

# def matriz_para_dicio(K, posicao):
#     dicio = {}
#     for linha in range(len(K)):
#         for coluna in range(len(K[linha])):
#             if posicao == 1:
#                 dicio[linha, coluna] = [(linha, coluna),K[linha][coluna]]
#             if posicao == 2:
#                 dicio[linha, coluna] = [(linha+2, coluna+2),K[linha][coluna]]
#             if posicao == 3:
#                 if linha == 0 or linha == 1:
#                     if coluna == 0 or coluna == 1:
#                         dicio[linha, coluna] = [(linha+4, coluna+4),K[linha][coluna]]
#                     elif coluna == 2 or coluna == 3:
#                         dicio[linha, coluna] = [(linha+4, coluna-2),K[linha][coluna]]
#                 elif linha == 2 or linha == 3:
#                     if coluna == 0 or coluna == 1:
#                         dicio[linha, coluna] = [(linha-2, coluna+4),K[linha][coluna]]
#                     elif coluna == 2 or coluna == 3:
#                         dicio[linha, coluna] = [(linha-2, coluna-2),K[linha][coluna]]
#     return dicio

# dicio1 = matriz_para_dicio(K1, 1)
# dicio2 = matriz_para_dicio(K2, 2)
# dicio3 = matriz_para_dicio(K3, 3)

def calculando_Kg(dicio1, dicio2, dicio3):
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
    return Kg


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
meuArquivo.write(f"Reacoes de apoio [N]\posicao{apoios}\posicao")
meuArquivo.write(f"\nDeslocamentos [m]\posicao{deslocamentos_expandido}\posicao")
meuArquivo.close()

