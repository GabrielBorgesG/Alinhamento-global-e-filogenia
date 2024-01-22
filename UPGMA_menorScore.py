import csv
from dataclasses import replace
import numpy as np

# Menor valor: localiza, dentro de uma matriz, o menor valor
def min_matriz(matriz):
    # Inicializa o menor valor como infinito (um nuumero muito grande)
    min_val = float("inf")
    x, y = -1, -1

    # Procura o menor valor
    for i in range(len(matriz)):
        for j in range(len(matriz[i])):
            if matriz[i][j] < min_val:
                min_val = matriz[i][j]
                x, y = i, j

    # Retorna as coordenadas do menor valor dentro da matriz
    return x, y

# Junta a identificacao das colunas das sequencias mais proximas
def Junta_Seq(nomes, a, b):
    # Organiza a ordem das sequencias pela identificacao
    if b < a:
        a, b = b, a

    # Junta as colunas na primeira coluna
    nomes[a] = "(" + nomes[a] + "," + nomes[b] + ")"

    # Remove a coluna extra
    del nomes[b]

# Junta colunas com as sequencias mais proximas
def Junta_coluna(matriz, a, b):
    # Organiza/ordena pelo indice
    if b < a:
        a, b = b, a

    # Recalcula e salva numa lista os valores da nova coluna
    linha = []
    for i in range(0, a):
        linha.append((matriz[a][i] + matriz[b][i])/2)
    matriz[a] = linha
    for i in range(a+1, b):
      matriz[i][a] = (matriz[i][a] + matriz[b][i])/2
    for i in range(b+1, len(matriz)):
        matriz[i][a] = (matriz[i][a] + matriz[i][b])/2
        # Remove a coluna redundante
        del matriz[i][b]

    # Remove a linha redundante
    del matriz[b]

# Inicializa os nomes de identificação das sequências
def ini_nomes(inicio, fim):
    nomes = []
    for i in range(ord(inicio), ord(fim)+1):
        nomes.append(chr(i))
    return nomes

# Integração das funções
def UPGMA(matriz_dist, nomes):
    while len(nomes) > 1:
        x, y = min_matriz(matriz_dist)
        Junta_coluna(matriz_dist, x, y)
        Junta_Seq(nomes, x, y)

    # Retorna o Newick
    return nomes[0]

# CÓDIGO

matriz_dist = [
[],
[-1303],
[-1580,	-1733],
[-1833,	-2139,	-1464],
[-1700,	-2125,	-1399,	-950],
[-1712,	-1913,	-1374,	-1247,	-1234],
[1254,	-1214,	-1607,	-1542,	-1987,	-1733],
[-732,	-1377,	-1636,	-1929,	-2153,	-1811,	-948],
[-1391,	-1534,	-1550,	-1589,	-1636,	-1502,	-1397,	-1411],
[-1470,	-1748,	-1444,	-1475,	-1497,	-1433,	-1483,	-1537,	-1504]
] # Obtido através do alinhamento no programa do Tabalho I

NOMES = ini_nomes("A", "J") # De A até H
print(NOMES)
newwick_final = UPGMA(matriz_dist, NOMES)

from Bio import Phylo
from Bio.Phylo.PhyloXML import Phylogeny
from io import StringIO

handle = StringIO(newwick_final)
tree = Phylo.read(handle, "newick")
Phylo.draw(tree)