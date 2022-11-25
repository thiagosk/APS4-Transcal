import numpy as np
from funcoesTermosol import importa
from numpy.linalg import det
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada_exemplo.xlsx')

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
        self.dicio = {}
    
    def get_posicao(self):
        return self.posicao

    def make_dicio(self, posicao_trelica):
        self.dicio = matriz_para_dicio(self.K, posicao_trelica)


def criacao_elementos():
    elementos = []  
    for i in range(nm):
        e = elemento(i+1, Inc[i][0], Inc[i][1], Inc[i][2], Inc[i][3])
        elementos.append(e)
    return elementos

elementos = criacao_elementos()


def criacao_trelicas(elementos):
    trelicas = []
    trelica_posicoes = []
    trelica_nos = []
    for element in elementos:
        trelica = []
        t_posicoes = []
        t_nos = []
        for e2 in elementos:
            if e2.no0 == element.nof:
                for e3 in elementos:
                    if e3.no0 == e2.nof or e3.nof == e2.nof:
                        if e3.nof == element.no0 or e3.no0 == element.no0:
                            trelica.append(element)
                            trelica.append(e2)
                            trelica.append(e3)
                            t_posicoes.append(element.posicao)
                            t_posicoes.append(e2.posicao)
                            t_posicoes.append(e3.posicao)
                            t_nos.append(int(element.no0))
                            t_nos.append(int(e2.no0))
                            if e3.no0 == element.no0:
                                t_nos.append(int(e3.nof))
                            else:
                                t_nos.append(int(e3.no0))
                            break
        if trelica != []:
            if sorted(t_posicoes) not in trelica_posicoes:
                trelicas.append(trelica)
                trelica_posicoes.append(t_posicoes)
                trelica_nos.append(t_nos)

    return trelicas, trelica_posicoes, trelica_nos

trelicas, trelica_posicoes, trelica_nos = criacao_trelicas(elementos)


def calcula_Kg(trelica, trelica_posicoes):

    # formando o dicio de cada elemento
    # for n_trelicas in range(len(trelicas)):
    #     for n_elementos in range(len(trelicas[n_trelicas])):
    #         for lista in trelica_posicoes:
    #             for n_posicoes in range(len(lista)):
    #                 if trelicas[n_trelicas][n_elementos].posicao == lista[n_posicoes]:
    #                     posicao_trelica = n_posicoes+1
    #         trelicas[n_trelicas][n_elementos].make_dicio(posicao_trelica)
    for i in range(len(trelica)):
        trelica[i].make_dicio(i+1)


    Kg = np.zeros((6,6))
    dicio1 = trelica[0].dicio
    dicio2 = trelica[1].dicio
    dicio3 = trelica[2].dicio
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

Kg_lista = []
for i in range(len(trelicas)):
    Kg_lista.append(calcula_Kg(trelicas[i], trelica_posicoes[i]))

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


def calcula_Pg(posicoes):
    Pg = []
    Pg_posicionados_indices = []
    for posicao in posicoes:
        if posicao != 1:
            indices = [2*(posicao-1), 2*(posicao-1)+1]
            if F[indices[0]] == 0:
                con = False
                for j in R:
                    if indices[0] == j[0]:
                        Pg.append("R")
                        con = True
                if con == False:
                    Pg.append(0)
            else:
                Pg.append(F[indices[0]][0])

            if F[indices[1]] == 0:
                c = False
                for j in R:
                    if indices[1] == j[0]:
                        Pg.append("R")
                        c = True
                if c == False:
                    Pg.append(0)
            else:
                Pg.append(F[indices[1]][0])
            Pg_posicionados_indices.append(indices[0])
            Pg_posicionados_indices.append(indices[1])
        else:
            if F[posicao-1] == 0:
                c = False
                for j in R:
                    if posicao-1 == j[0]:
                        Pg.append("R")
                        c = True
                if c == False:
                    Pg.append(0)
            else:
                Pg.append(F[posicao-1][0])

            if F[posicao][0] == 0:
                con = False
                for j in R:
                    if posicao == j[0]:
                        Pg.append("R")
                        con = True
                if con == False:
                    Pg.append(0)
            else:
                Pg.append(F[posicao][0])
            Pg_posicionados_indices.append(posicao-1)
            Pg_posicionados_indices.append(posicao)
    return Pg, Pg_posicionados_indices


def filtra_Pg(Pg, Pg_posicionados_indices):
    lista_indices_nao_apoio_Pg = []
    lista_indices_apoio_Pg = []
    Pg_filtrado = []
    l_nao_apoio_posicoes_indice = []
    for i in range(len(Pg)):
        if Pg[i] != "R":
            lista_indices_nao_apoio_Pg.append(i)
            Pg_filtrado.append(Pg[i])
            l_nao_apoio_posicoes_indice.append(Pg_posicionados_indices[i])
        else:
            lista_indices_apoio_Pg.append(i)
    return Pg_filtrado, lista_indices_apoio_Pg, lista_indices_nao_apoio_Pg, l_nao_apoio_posicoes_indice
    

Pg_lista = []
Pg_filtrado_lista = []
lista_indices_apoio_Pg_lista = []
lista_indices_nao_apoio_Pg_lista = []
lista_l_nao_apoio_posicoes_indice = []
for posicoes in trelica_nos:
    Pg, Pg_posicionados_indices = calcula_Pg(posicoes)
    Pg_filtrado, lista_indices_apoio_Pg, lista_indices_nao_apoio_Pg, l_nao_apoio_posicoes_indice = filtra_Pg(Pg, Pg_posicionados_indices)
    Pg_lista.append(Pg)
    Pg_filtrado_lista.append(Pg_filtrado)
    lista_indices_apoio_Pg_lista.append(lista_indices_apoio_Pg)
    lista_indices_nao_apoio_Pg_lista.append(lista_indices_nao_apoio_Pg)
    lista_l_nao_apoio_posicoes_indice.append(l_nao_apoio_posicoes_indice)


def filtra_Kg(Kg, lista_indices_nao_apoio_Pg, Pg, Pg_total):
    posicao_em_Pg_total = []
    for linha in range(len(Kg)):
        if linha in lista_indices_nao_apoio_Pg:
            posicao_em_Pg_total.append(linha)

    Kg_filtrado_linha = []
    posicao_em_Pg = []
    for linha in range(len(Kg)):
        if linha in lista_indices_nao_apoio_Pg:
            Kg_filtrado_linha.append(Kg[linha])
            posicao_em_Pg.append(linha)

    Pg_filtrado_trelica = []
    for posicao in posicao_em_Pg:
        Pg_filtrado_trelica.append(Pg[posicao])

    Kg_filtrado = []
    for linha in range(len(Kg_filtrado_linha)):
        l = []
        for coluna in range(len(Kg_filtrado_linha[linha])):
            if coluna in lista_indices_nao_apoio_Pg:
                l.append(Kg_filtrado_linha[linha][coluna])
        Kg_filtrado.append(l)
    return Kg_filtrado, posicao_em_Pg, Pg_filtrado_trelica

Kg_filtrado_lista = []
posicao_em_Pg_lista = []
Pg_filtrado_trelica_lista = []
for i in range(len(Kg_lista)):
    Kg_filtrado, posicao_em_Pg, Pg_filtrado_trelica = filtra_Kg(Kg_lista[i], lista_indices_nao_apoio_Pg_lista[i], Pg_lista[i], Pg_total)
    Kg_filtrado_lista.append(Kg_filtrado)
    posicao_em_Pg_lista.append(posicao_em_Pg)
    Pg_filtrado_trelica_lista.append(Pg_filtrado_trelica)

def calculo_deslocamentos(Kg_filtrado, Pg_filtrado_trelica):
    multiplicador_das_incognitas = []
    apoios = Pg_filtrado_trelica

    if det(Kg_filtrado) == 0:
        return np.zeros((len(Pg_filtrado_trelica), len(Pg_filtrado_trelica)))

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
    return deslocamentos


deslocamentos_dicio = {}
for j in range(len(Kg_filtrado_lista)):
    deslocamentos = calculo_deslocamentos(Kg_filtrado_lista[j], Pg_filtrado_trelica_lista[j])
    for i in range(len(lista_l_nao_apoio_posicoes_indice[j])):
        deslocamentos_dicio[lista_l_nao_apoio_posicoes_indice[j][i]] = deslocamentos[i]


deslocamentos_expandido = np.zeros((len(Pg_total),1))
contagem = np.arange(0, len(Pg_total))
for i in deslocamentos_dicio:
    try:
        if deslocamentos_dicio[i]:
            deslocamentos_expandido[i] = deslocamentos_dicio[i]
        else:
            deslocamentos_expandido[i] = 0
    except:
        deslocamentos_expandido[i] = 0
    # if deslocamentos_dicio[i].is_array():
    #     deslocamentos_expandido[i] = 0
    # else:
    #     if deslocamentos_dicio[i]:
    #         deslocamentos_expandido[i] = deslocamentos_dicio[i]
    #     else:
    #         deslocamentos_expandido[i] = 0

# # Obtendo as reações de apoio estrutural

apoios =  []
for i in range(len(Kg_lista)):
    for linha in range(len(Kg_lista[i])):
        calculo = [0]
        if linha in lista_indices_apoio_Pg_lista[i]:
            for coluna in range(len(Kg_lista[i][linha])):
                calculo += Kg_lista[i][linha][coluna] * deslocamentos_expandido[coluna]
            apoios.append([round(calculo[0],1)])


# # Escrevendo o arquivo de saída

meuArquivo = open('saida.txt', 'w')
meuArquivo.write(f"Reacoes de apoio [N]\n{apoios}\n")
meuArquivo.write(f"\nDeslocamentos [m]\n{deslocamentos_expandido}\n")
meuArquivo.close()

