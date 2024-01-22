#cd '.\Desktop\INFORMÁTICA\BioComp\Trab II\'

import csv
import numpy as np
import pandas as pd

#pd.set_option('max_columns', None)
matriz_neighbor_joining = np.zeros((8, 8)) # A-H X A-H

print(matriz_neighbor_joining)
# Preenchimento da matriz

matriz_neighbor_joining[0][1] = 2.2266620
matriz_neighbor_joining[0][2] = 3.3353841	
matriz_neighbor_joining[0][3] = 3.0926941
matriz_neighbor_joining[0][4] = 1.8896296
matriz_neighbor_joining[0][5] = 4.9589915	
matriz_neighbor_joining[0][6] = 4.3943885	
matriz_neighbor_joining[0][7] = 5.4122139
matriz_neighbor_joining[1][0] = matriz_neighbor_joining[0][1]
matriz_neighbor_joining[1][2] = 1.1567076
matriz_neighbor_joining[1][3] = 0.9140176
matriz_neighbor_joining[1][4] = 1.3807486
matriz_neighbor_joining[1][5] = 4.4501105
matriz_neighbor_joining[1][6] = 3.8855074
matriz_neighbor_joining[1][7] = 4.9033329
matriz_neighbor_joining[2][0] = matriz_neighbor_joining[0][2]
matriz_neighbor_joining[2][1] = matriz_neighbor_joining[1][2]
matriz_neighbor_joining[2][3] = 1.8198296	
matriz_neighbor_joining[2][4] = 2.4894707
matriz_neighbor_joining[2][5] = 5.5588327
matriz_neighbor_joining[2][6] = 4.9942296
matriz_neighbor_joining[2][7] = 6.0120551
matriz_neighbor_joining[3][0] = matriz_neighbor_joining[0][3]
matriz_neighbor_joining[3][1] = matriz_neighbor_joining[1][3]
matriz_neighbor_joining[3][2] = matriz_neighbor_joining[2][3]
matriz_neighbor_joining[3][4] = 2.2467807	
matriz_neighbor_joining[3][5] = 5.3161427
matriz_neighbor_joining[3][6] = 4.7515396
matriz_neighbor_joining[3][7] = 5.7693651
matriz_neighbor_joining[4][0] = matriz_neighbor_joining[0][4]
matriz_neighbor_joining[4][1] = matriz_neighbor_joining[1][4]
matriz_neighbor_joining[4][2] = matriz_neighbor_joining[2][4]
matriz_neighbor_joining[4][3] = matriz_neighbor_joining[3][4]
matriz_neighbor_joining[4][5] = 3.3471233
matriz_neighbor_joining[4][6] = 2.7825202
matriz_neighbor_joining[4][7] = 3.8003457
matriz_neighbor_joining[5][0] = matriz_neighbor_joining[0][5]
matriz_neighbor_joining[5][1] = matriz_neighbor_joining[1][5]
matriz_neighbor_joining[5][2] = matriz_neighbor_joining[2][5]
matriz_neighbor_joining[5][3] = matriz_neighbor_joining[3][5]
matriz_neighbor_joining[5][4] = matriz_neighbor_joining[4][5]
matriz_neighbor_joining[5][6] = 1.1966232
matriz_neighbor_joining[5][7] = 2.2144487
matriz_neighbor_joining[6][0] = matriz_neighbor_joining[0][6]
matriz_neighbor_joining[6][1] = matriz_neighbor_joining[1][6]
matriz_neighbor_joining[6][2] = matriz_neighbor_joining[2][6]
matriz_neighbor_joining[6][3] = matriz_neighbor_joining[3][6]
matriz_neighbor_joining[6][4] = matriz_neighbor_joining[4][6]
matriz_neighbor_joining[6][5] = matriz_neighbor_joining[5][6]
matriz_neighbor_joining[6][7] = 1.2176200
matriz_neighbor_joining[7][0] = matriz_neighbor_joining[0][7]
matriz_neighbor_joining[7][1] = matriz_neighbor_joining[1][7]
matriz_neighbor_joining[7][2] = matriz_neighbor_joining[2][7]
matriz_neighbor_joining[7][3] = matriz_neighbor_joining[3][7]
matriz_neighbor_joining[7][4] = matriz_neighbor_joining[4][7]
matriz_neighbor_joining[7][5] = matriz_neighbor_joining[5][7]
matriz_neighbor_joining[7][6] = matriz_neighbor_joining[6][7]

# Cálculo da coluna de Distância total
vetor_dist_total = []
for x in range (8):
    n = 0
    for y in range (8):
        n = n + matriz_neighbor_joining[x][y]
    vetor_dist_total.append(n)

matriz_distancia = np.zeros((8, 8)) # A-H X A-H

#M(A,B) = d(A, B) - (r(A) + r(B))/(n taxons - 2)

#matriz_distancia[0][0] = matriz_neighbor_joining[0][1] - (vetor_dist_total[0] + vetor_dist_total[1])/(8-2)
#print(matriz_distancia[0][0]) # AB 2,2266620 - (25,3099637 + 18,9170866)/6 = -5.14451305

#matriz_distancia[3][3] = matriz_neighbor_joining[3][4] - (vetor_dist_total[3] + vetor_dist_total[4])/(8-2)
#print(matriz_distancia[3][3]) # DE 2,2467807 - (23,9103694 + 17,9366188)/6 = -4,72771733

#matriz_distancia[6][6] = matriz_neighbor_joining[6][7] - (vetor_dist_total[6] + vetor_dist_total[7])/(8-2)
#print(matriz_distancia[6][6]) # GH 1,2176200 - (23,2224285 + 29,3293814)/6 = -7,54101498

for k in range(8): # linhas
    for w in range(8): # colunas
        if(k != w):
            matriz_distancia[k][w] = matriz_neighbor_joining[k][w] - (vetor_dist_total[k] + vetor_dist_total[w])/(8-2)

# PRIMEIRA RODADA
# Encontra a posição de menor distância
menor = 0
coluna, linha = 0, 0
for k in range(8): # linhas
    for w in range(8): # colunas
        if(matriz_distancia[k][w] < menor):
            menor = matriz_distancia[k][w]
            coluna = w # Guarda coluna
            linha = k # Guarda linha

S1 = (matriz_distancia[linha][coluna]/2) + ((vetor_dist_total[linha] - vetor_dist_total[coluna])/(2*6))
S2 = matriz_distancia[linha][coluna] - S1
# S1 = -7,54101498/2 + (29,3293814 - 23,2224285)/(2*(8-2)) = -4,279420231
# S2 = -3,26159474

# Calcular as novas distâncias:
#d(S, outra seq) = [d(Seq1, outra seq) + d(Seq2, outra seq) - d(Seq1, Seq2)]/2

# Transformando a matriz para trabalhar com Pandas
#matriz_distancia = np.append(matriz_distancia, [[1], [1], [1], [1], [1], [1], [1]], axis=2)
MD = pd.DataFrame(matriz_distancia)
MD = MD.rename(columns={0:'A',1:'B',2:'C',3:'D',4:'E',5:'F',6:'G', 7: 'H'}) # Muda os nomes das colunas
MD.index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', ] # Muda os nomes das linhas

MD.insert(8,"GH", [1.0, 1, 1, 1, 1, 1, 1, 1], True) # Adicionando coluna
matriz_distancia = np.append(matriz_distancia, [[1], [1], [1], [1], [1], [1], [1], [1]], axis=1) # Adicionando coluna
#matriz_distancia = np.append(matriz_distancia, [[1], [1], [1], [1], [1], [1], [1], [1], [1]]) # Adiciona linha
#MD.loc[len(MD.index )] = [1, 1, 1, 1, 1, 1, 1, 1, 1] # Adiciona linha

#nova_linha = {'GH':1, 'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E':1, 'F':1, 'G':1}
#MD = MD.append(nova_linha, ignore_index=True, sort = False) # Adicionando linha
#MD.loc['GH'] = [1, 1, 1, 1, 1, 1, 1, 1] # Adicionando linha
print(matriz_distancia)
# Lista que será coluna:
coluna_nova = []
newick = ""
newick2 = ""
newick = newick + "(G, H)"

for i in range(8): # Existem 8 linhas por enquanto
    coluna_nova.append([matriz_distancia[i][coluna] + matriz_distancia[i][linha] - matriz_neighbor_joining[linha][coluna]][0]/2.0) # Salva os valores da nova coluna

for i in range(8):
    MD["GH"][i] = coluna_nova[i] #Atribui os novos valores
    matriz_distancia[i][8] = coluna_nova[i] #Atribui os novos valores
    
MD = MD.drop(MD.index[linha])
MD = MD.drop(MD.index[coluna-1])
MD = MD.drop(columns=["G"])
MD = MD.drop(columns=["H"])

matriz_distancia = np.delete(matriz_distancia, coluna, 0) # Deleta a linha H
matriz_distancia = np.delete(matriz_distancia, linha, 0) # Deleta a linha G
matriz_distancia = np.delete(matriz_distancia, coluna, 1) # Deleta a coluna H
matriz_distancia = np.delete(matriz_distancia, linha, 1) # Deleta a coluna G
print(MD)
# SEGUNDA RODADA
menor = 0
coluna, linha = 0, 0
for k in range(6): # linhas
    for w in range(7): # colunas
        if(matriz_distancia[k][w] < menor):
            menor = matriz_distancia[k][w]
            coluna = w # Guarda coluna
            linha = k # Guarda linha
print(linha, coluna, matriz_distancia[linha][coluna])
# Encontrou que a menor distância é entre GH e F

S1 = (matriz_distancia[linha][coluna]/2) + ((vetor_dist_total[linha] - vetor_dist_total[coluna])/(2*6))
S2 = matriz_distancia[linha][coluna] - S1

MD.insert(7,"FGH", [1.0, 1, 1, 1, 1, 1], True) # Adicionando coluna
matriz_distancia = np.append(matriz_distancia, [[1], [1], [1], [1], [1], [1]], axis=1) # Adicionando coluna
newick = "(" + newick + ", F)"

coluna_nova = []
for i in range(6):
    coluna_nova.append([matriz_distancia[i][coluna] + matriz_distancia[i][linha] - matriz_neighbor_joining[linha][coluna]][0]/2.0)

for i in range(6):
    MD["FGH"][i] = coluna_nova[i]
    matriz_distancia[i][7] = coluna_nova[i] #Atribui os novos valores

MD = MD.drop(MD.index[linha])
MD = MD.drop(columns=["F"])
MD = MD.drop(columns=["GH"])

matriz_distancia = np.delete(matriz_distancia, linha, 0)
matriz_distancia = np.delete(matriz_distancia, coluna, 1)
matriz_distancia = np.delete(matriz_distancia, linha, 1)
print(MD)
# TERCEIRA RODADA
menor = 0
coluna, linha = 0, 0
for k in range(5): # linhas
    for w in range(6): # colunas
        if(matriz_distancia[k][w] < menor):
            menor = matriz_distancia[k][w]
            coluna = w # Guarda coluna
            linha = k # Guarda linha
print(linha, coluna, matriz_distancia[linha][coluna])
# Encontrou que a menor distância é entre GH e F

S1 = (matriz_distancia[linha][coluna]/2) + ((vetor_dist_total[linha] - vetor_dist_total[coluna])/(2*6))
S2 = matriz_distancia[linha][coluna] - S1

MD.insert(6,"CD", [1.0, 1, 1, 1, 1], True) # Adicionando coluna
matriz_distancia = np.append(matriz_distancia, [[1], [1], [1], [1], [1]], axis=1) # Adicionando coluna
newick2 = newick2 + "(C, D)"

coluna_nova = []
for i in range(5):
    coluna_nova.append([matriz_distancia[i][coluna] + matriz_distancia[i][linha] - matriz_neighbor_joining[linha][coluna]][0]/2.0)
print(coluna_nova)

for i in range(5):
    MD["CD"][i] = coluna_nova[i]
    matriz_distancia[i][6] = coluna_nova[i] # Atribui os novos valores
print(MD)
MD = MD.drop(MD.index[coluna], axis = 0)
MD = MD.drop(MD.index[linha], axis = 0)
MD = MD.drop(columns=["C"], axis = 1)
MD = MD.drop(columns=["D"], axis = 1)

matriz_distancia = np.delete(matriz_distancia, coluna, 0)
matriz_distancia = np.delete(matriz_distancia, linha, 0)
matriz_distancia = np.delete(matriz_distancia, coluna, 1)
matriz_distancia = np.delete(matriz_distancia, linha, 1)

print(MD)
print(matriz_distancia)

# QUARTA RODADA
menor = 0
coluna, linha = 0, 0
for k in range(3): # linhas
    for w in range(5): # colunas
        if(matriz_distancia[k][w] < menor):
            menor = matriz_distancia[k][w]
            coluna = w # Guarda coluna
            linha = k # Guarda linha
print(linha, coluna, matriz_distancia[linha][coluna])
# Encontrou que a menor distância é entre GH e F
S1 = (matriz_distancia[linha][coluna]/2) + ((vetor_dist_total[linha] - vetor_dist_total[coluna])/(2*6))
S2 = matriz_distancia[linha][coluna] - S1

MD.insert(5,"BCD", [1.0, 1, 1], True) # Adicionando coluna
matriz_distancia = np.append(matriz_distancia, [[1], [1], [1]], axis=1) # Adicionando coluna
newick2 = "(B, " + newick2 + ")"

coluna_nova = []
for i in range(3):
    coluna_nova.append([matriz_distancia[i][coluna] + matriz_distancia[i][linha] - matriz_neighbor_joining[linha][coluna]][0]/2.0)
print(coluna_nova)

for i in range(3):
    MD["BCD"][i] = coluna_nova[i]
    matriz_distancia[i][5] = coluna_nova[i] # Atribui os novos valores

MD = MD.drop(MD.index[linha])
MD = MD.drop(columns=["CD"])
MD = MD.drop(columns=["B"])

matriz_distancia = np.delete(matriz_distancia, linha, 0)
matriz_distancia = np.delete(matriz_distancia, coluna, 1)
matriz_distancia = np.delete(matriz_distancia, linha, 1)

# QUINTA RODADA
menor = 0
coluna, linha = 0, 0
for k in range(2): # linhas
    for w in range(4): # colunas
        if(matriz_distancia[k][w] < menor):
            menor = matriz_distancia[k][w]
            coluna = w # Guarda coluna
            linha = k # Guarda linha
print(linha, coluna, matriz_distancia[linha][coluna])
# Encontrou que a menor distância é entre GH e F
S1 = (matriz_distancia[linha][coluna]/2) + ((vetor_dist_total[linha] - vetor_dist_total[coluna])/(2*6))
S2 = matriz_distancia[linha][coluna] - S1

MD.insert(4,"ABCD", [1.0, 1], True) # Adicionando coluna
matriz_distancia = np.append(matriz_distancia, [[1], [1]], axis=1) # Adicionando coluna
newick2 = "(A, " + newick2 + ")"

coluna_nova = []
for i in range(2):
    coluna_nova.append([matriz_distancia[i][coluna] + matriz_distancia[i][linha] - matriz_neighbor_joining[linha][coluna]][0]/2.0)
print(coluna_nova)

for i in range(2):
    MD["ABCD"][i] = coluna_nova[i]
    matriz_distancia[i][4] = coluna_nova[i] # Atribui os novos valores

MD = MD.drop(MD.index[linha])
MD = MD.drop(columns=["BCD"])
MD = MD.drop(columns=["A"])

matriz_distancia = np.delete(matriz_distancia, linha, 0)
matriz_distancia = np.delete(matriz_distancia, coluna, 1)
matriz_distancia = np.delete(matriz_distancia, linha, 1)

# SEXTA RODADA
menor = 0
coluna, linha = 0, 0
for k in range(1): # linhas
    for w in range(3): # colunas
        if(matriz_distancia[k][w] < menor):
            menor = matriz_distancia[k][w]
            coluna = w # Guarda coluna
            linha = k # Guarda linha
print(linha, coluna, matriz_distancia[linha][coluna])
# Encontrou que a menor distância é entre GH e F
S1 = (matriz_distancia[linha][coluna]/2) + ((vetor_dist_total[linha] - vetor_dist_total[coluna])/(2*6))
S2 = matriz_distancia[linha][coluna] - S1

MD.insert(3,"ABCDE", [1.0], True) # Adicionando coluna
matriz_distancia = np.append(matriz_distancia, [[1]], axis=1) # Adicionando coluna
newick2 = "(" + newick2 + ", E)"

coluna_nova = []
for i in range(1):
    coluna_nova.append([matriz_distancia[i][coluna] + matriz_distancia[i][linha] - matriz_neighbor_joining[linha][coluna]][0]/2.0)
print(coluna_nova)

for i in range(1):
    MD["ABCDE"][i] = coluna_nova[i]
    matriz_distancia[i][3] = coluna_nova[i] # Atribui os novos valores

MD = MD.drop(columns=["ABCD"])
MD = MD.drop(columns=["E"])

matriz_distancia = np.delete(matriz_distancia, coluna, 1)
matriz_distancia = np.delete(matriz_distancia, linha, 1)

# SÉTIMA RODADA

newick = "(" + newick2 + newick + ")"

print(MD)
print(newick)
print(matriz_distancia)

from Bio import Phylo
from Bio.Phylo.PhyloXML import Phylogeny
from io import StringIO

handle = StringIO(newick)
tree = Phylo.read(handle, "newick")
Phylo.draw(tree)