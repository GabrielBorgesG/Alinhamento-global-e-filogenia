import numpy as np
#cd '.\Desktop\INFORMÁTICA\BioComp\Trab II\dataset_1\'

# CONSTANTES
MATCH = 1
MISMATCH = -1
GAP = -2

# VARIÁVEIS DE CÁLCULO
Cal1 = 0 # D(i-1)(j-1) match
Cal2 = 0 # D(i-1)(j-1) mismatch
Cal3 = 0 # D(i-1)(j) gap
Cal4 = 0 # D(i)(j-1) gap

# CÓDIGO
print("Alinhamento global de sequências\n")

# Entrada das sequências
with open('seq_9.fasta.fna', 'r') as f:
    seq1 = f.read()
with open('seq_10.fasta.fna', 'r') as f:
    seq2 = f.read()

# Cálculo para as dimensões da matriz
tam_seq1 = len(seq1) + 1 # Colunas
tam_seq2 = len(seq2) + 1 # Linhas

# Criação da matriz iniciada com zeros e preenchimento da primeira linha e da primeira coluna com valores de gaps
matriz_principal = np.zeros((tam_seq2, tam_seq1))
for k in range (1, tam_seq2):
    matriz_principal[k][0] = int(matriz_principal[k-1][0] + GAP)
for w in range (1, tam_seq1):
    matriz_principal[0][w] = int(matriz_principal[0][w-1] + GAP)

# Calcula a matriz
for x in range (1, tam_seq2):
    for y in range (1, tam_seq1):
        Cal1 = matriz_principal[x-1][y-1] + MATCH
        Cal2 = matriz_principal[x-1][y-1] + MISMATCH
        Cal3 = matriz_principal[x-1][y] + GAP
        Cal4 = matriz_principal[x][y-1] + GAP
        # Confere se tem match
        if(seq1[y-1] == seq2[x-1]): # Se tiver match, Cal1 é válido
            matriz_principal[x][y] = int(max(Cal1, Cal3, Cal4)) # Encontra o maior valor e armazena nesta posição da matriz
        else: # Se não, Cal1 é inválido
            matriz_principal[x][y] = int(max(Cal2, Cal3, Cal4)) # Encontra o maior valor e armazena nesta posição da matriz

print("\n")
print("O score do alinhamento é: ", matriz_principal[tam_seq2-1][tam_seq1-1])
print("\n")