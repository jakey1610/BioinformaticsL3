import math
from itertools import product
import numpy as np
from sympy import pprint, Matrix
#When referring to X or Y the index should be one lower than for tables.
def scoreFunc(x,y):
    score = 0
    if len(x)!=len(y):
        return "Length Error"
    for i in range(len(x)):
        if x[i] == y[i]:
            score += 1
        elif x[i] == "-" or y[i] == "-":
            score -= 2
        else:
            score -= 1
    return score

def reconstSeq(d,s,X,Y,pX,pY):
    ind1 = []
    ind2 = []
    i,j,sc = pX,pY,0
    for row in range(len(s)):
        for score in range(len(s[row])):
            if s[row][score] > sc:
                sc = s[row][score]
                i = row
                j = score

    while d[i][j]!="E":
        if d[i][j] == "D":
            i-=1
            j-=1
            if X[i] == Y[j]:
                ind1.append(i)
                ind2.append(j)
        elif d[i][j] == "L":
            i-=1
        elif d[i][j] == "U":
            j-=1

    return ind1, ind2, sc

#### DYNAMIC PROGRAMMING ####
def dynprog(alphabet, scoring_matrix, sequence1, sequence2):
    pSq1 = len(sequence1)
    pSq2 = len(sequence2)
    scoreMat = [[-1 for sequence1 in range(pSq2+1)] for sequence2 in range(pSq1+1)]
    dirMat = [["U" for sequence1 in range(pSq2+1)] for sequence2 in range(pSq1+1)]
    dirMat[0][0] = "E"
    scoreMat[0][0] = 0

    for i in range(1,pSq1+1):
        indexChar = alphabet.index(sequence1[i-1])
        scoreMat[i][0] = max(0,scoreMat[i-1][0] + scoring_matrix[indexChar][-1])
        dirMat[i][0] = "E"
    for i in range(1,pSq2+1):
        indexChar = alphabet.index(sequence2[i-1])
        scoreMat[0][i] = max(0,scoreMat[0][i-1] + scoring_matrix[-1][indexChar])
        dirMat[0][i] = "E"

    for j in range(pSq2):
        for i in range(pSq1):
            indexChars = [alphabet.index(sequence1[i]), alphabet.index(sequence2[j])]
            scoreMat[i+1][j+1] = max(scoreMat[i][j] + scoring_matrix[indexChars[0]][indexChars[1]], scoreMat[i][j+1] + scoring_matrix[indexChars[0]][-1], scoreMat[i+1][j] + scoring_matrix[-1][indexChars[1]], 0)
            if scoreMat[i+1][j+1] == scoreMat[i][j] + scoring_matrix[indexChars[0]][indexChars[1]]:
                dirMat[i+1][j+1] = "D"
            elif scoreMat[i+1][j+1] == scoreMat[i][j+1] + scoring_matrix[indexChars[0]][-1]:
                dirMat[i+1][j+1] = "L"
            elif scoreMat[i+1][j+1] == scoreMat[i+1][j] + scoring_matrix[-1][indexChars[1]]:
                dirMat[i+1][j+1] = "U"
            elif scoreMat[i+1][j+1] == 0:
                dirMat[i+1][j+1] = "E"
    ind1,ind2,score = reconstSeq(dirMat,scoreMat,sequence1,sequence2,pSq1,pSq2)
    return score, ind1[::-1], ind2[::-1]

#### DYNAMIC PROGRAMMING IN LINEAR SPACE ####
def NWScore(alphabet, scoring_matrix, sequence1, sequence2):
    pSq1 = len(sequence1)
    pSq2 = len(sequence2)
    F = [[0 for x in range(2)] for y in range(pSq1+1)]
    B = [[0 for x in range(2)] for y in range(pSq1+1)]
    #Initialise F array
    for i in range(1,pSq1+1):
        indexChar = alphabet.index(sequence1[i-1])
        F[i][0] = max(F[i-1][0] + scoring_matrix[indexChar][-1],0)
    #Populate F array
    for j in range(pSq2):
        for i in range(pSq1):
            indexChars = [alphabet.index(sequence1[i]), alphabet.index(sequence2[j])]
            F[i+1][(j+1)%2] = max(F[i][(j)%2] + scoring_matrix[indexChars[0]][indexChars[1]], F[i][(j+1)%2] + scoring_matrix[indexChars[0]][-1], F[i+1][j%2] + scoring_matrix[-1][indexChars[1]], 0)

    #Initialise B array
    for j in range(0,pSq1+1):
        indexChars = [alphabet.index(sequence1[j-1]), alphabet.index(sequence2[pSq2-1])]
        B[j][1] = max(scoring_matrix[indexChars[0]][indexChars[1]],0)
    #Populate B array
    for j in range(pSq2,0,-1):
        for i in range(pSq1,0,-1):
            indexChars = [alphabet.index(sequence1[i-1]), alphabet.index(sequence2[j-1])]
            B[i-1][(j-1)%2] = max(B[i][j%2] + scoring_matrix[indexChars[0]][indexChars[1]], B[i-1][j%2] + scoring_matrix[indexChars[0]][-1], B[i][(j-1)%2] + scoring_matrix[-1][indexChars[1]],0)
    return F,B

def dynproglin(alphabet, scoring_matrix, seq1,seq2,offsetx=0,offsety=0):
    pSq1 = len(seq1)
    pSq2 = len(seq2)
    if pSq1 >= pSq2:
        sequence1 = seq1
        sequence2 = seq2
    else:
        sequence1 = seq2
        sequence2 = seq1
        temp = pSq1
        pSq1 = pSq2
        pSq2 = temp
    ind1 = []
    ind2 = []
    if pSq2==1:
        l1 = []
        l2 = []
        ind = 0
        for i in range(0,pSq1):
            ind = i
            if sequence1[ind] == sequence2[0]:
                l1 = [ind+offsetx]
                l2 = [offsety]
                break

        return (0, l1,l2)
    elif pSq1 == 1:
        l1 = []
        l2 = []
        ind = 0
        for i in range(0,pSq2):
            ind = i
            if sequence1[0] == sequence2[ind]:
                l1 = [offsetx]
                l2 = [ind+offsety]
                break
        return (0, l1, l2)
    elif pSq1 == 0 or pSq2 == 0:
        return (0,[],[])
    else:
        F,B = NWScore(alphabet, scoring_matrix, sequence1, sequence2)

        fH, b = NWScore(alphabet, scoring_matrix, sequence2, sequence1[:math.floor(pSq1/2)])
        f, bH = NWScore(alphabet, scoring_matrix, sequence2, sequence1[math.floor(pSq1/2):])
        pos = 0
        inter = -1
        for i in range(min(len(fH),len(bH))):
            inter = max(inter, fH[i][1] + bH[i][0])
            if inter == fH[i][1] + bH[i][0]:
                pos = i
        I = dynproglin(alphabet, scoring_matrix, sequence1[:math.floor(pSq1/2)],sequence2[:pos], offsetx=offsetx, offsety=offsety)
        J = dynproglin(alphabet, scoring_matrix, sequence1[math.floor(pSq1/2):],sequence2[pos:],offsetx=offsetx+math.floor(pSq1/2), offsety=offsety+pos)
        ind1 += I[1] + J[1]
        ind2 += I[2] + J[2]
        sc = 0
        for row in range(len(F)):
            for col in range(len(F[row])):
                if F[row][col] > sc:
                    sc = F[row][col]

        return sc,ind1,ind2

#### HEURISTIC ALGORITHM ####
def heuralign(alphabet, scoring_matrix, seq1, seq2):
    #sequence1 is the shorter of the two sequences
    switch = False

    sequence1 = seq1
    sequence2 = seq2

    pSq1 = len(sequence1)
    pSq2 = len(sequence2)
    lnth = 2
    k = 32
    index_table = dict()
    l=[]
    for i in range(lnth,max(0,lnth-4),-1):
        l += list(product(alphabet, repeat=i))
    l = [''.join(x) for x in l]
    for seq in l:
        index_table[seq] = []

    #populate the index table with small sequence indices
    for i in range(lnth,max(0,lnth-4),-1):
        ind = 0
        while ind < len(sequence1)-(i-1):
            if sequence1[ind:ind+i] in index_table.keys():
                index_table[sequence1[ind:ind+i]].append(ind)
            ind+=1

    #march on the query sequence2 checking for matches
    pot_matches = []
    for i in range(lnth,max(0,lnth-4),-1):
        ind = 0
        while ind < len(sequence2)-(i-1):
            if index_table[sequence2[ind:ind+i]] == []:
                pass
            else:
                for match in index_table[sequence2[ind:ind+i]]:
                    pot_matches.append((ind,match))
            ind+=1

    #This gets only one match from each diagonal
    pM = []
    pD = []
    for i in range(len(pot_matches)):
        if abs(pot_matches[i][0]-pot_matches[i][1]) not in pD:
            pM.append(pot_matches[i])
            pD.append(abs(pot_matches[i][0]-pot_matches[i][1]))

    #0th element is sequence2 and 1st element is sequence1
    pot_matches = pM
    score_matches = []
    #get the 7 best diagonals
    for p in pot_matches:
        score_matches.append(scoring_matrix[alphabet.index(sequence1[p[1]])][alphabet.index(sequence2[p[0]])])
    pot_matches = [x for _,x in sorted(zip(score_matches,pot_matches))][:7]

    #pot_diagonals will be a list of the scores on each diagonal
    #pot_ind will be the list of indices of matches
    pot_diagonals = []
    pot_ind = []
    for match in pot_matches:
        bMatch = bandedDP(alphabet, scoring_matrix, sequence1, sequence2, match, k)
        pot_diagonals.append(bMatch[0])
        pot_ind.append((bMatch[1],bMatch[2]))
    #return the gapped diagonal which has the best score
    score = np.amax(pot_diagonals)
    result = np.where(pot_diagonals == np.amax(pot_diagonals))
    result = list(np.array(result).astype('int16'))

    if switch:
        pot_ind[result[0][0]] = pot_ind[result[0][0]][::-1]
    return score, pot_ind[result[0][0]]



def bandedDP(alphabet, scoring_matrix, sequence1, sequence2, m, k):
    pSq1 = len(sequence1)
    pSq2 = len(sequence2)
    # scoreMat = list(np.zeros((pSq1+1,pSq2+1)).astype('int16'))
    scoreMat = dict()
    maxScore = -1
    maxInd = [-1,-1]
    diagonal = m[0]-m[1]
    if diagonal >= 0:
        i = diagonal - k
        j = 0
    else:
        i = 0
        j = -1*(diagonal) - k
    while i <= pSq1 + k and j <= pSq2 + k:
        for x in range(j-k,j+k):
            if x >= 0 and x <= pSq2 and i >= 0 and i<= pSq1:
                indexChars = [alphabet.index(sequence1[i-1]), alphabet.index(sequence2[x-1])]
                if x != 0 and i != 0:
                    # dS = scoreMat[i-1][x-1] + scoring_matrix[indexChars[0]][indexChars[1]]
                    try:
                        dS = scoreMat[(i-1,x-1)] + scoring_matrix[indexChars[0]][indexChars[1]]
                    except:
                        dS = scoring_matrix[indexChars[0]][indexChars[1]]
                else:
                    dS = -1
                if x != j-k:
                    # lS = scoreMat[i][x-1] + scoring_matrix[-1][indexChars[1]]
                    try:
                        lS = scoreMat[(i,x-1)] + scoring_matrix[-1][indexChars[1]]
                    except:
                        lS = scoring_matrix[-1][indexChars[1]]
                else:
                    lS = -1
                if x != j+k:
                    # uS = scoreMat[i-1][x] + scoring_matrix[indexChars[0]][-1]
                    try:
                        uS = scoreMat[(i-1,x)] + scoring_matrix[indexChars[0]][-1]
                    except:
                        uS = scoring_matrix[indexChars[0]][-1]
                else:
                    uS = -1
                scoreMat[(i,x)] = max(dS,lS,uS,0)
                if scoreMat[(i,x)] > maxScore:
                    maxScore = scoreMat[(i,x)]
                    maxInd = (i,x)
        i+=1
        j+=1
    i = maxInd[0]
    j = maxInd[1]
    ind1 = []
    ind2 = []

    while True:
        indexChars = [alphabet.index(sequence1[i-1]), alphabet.index(sequence2[j-1])]
        try:
            scoreMat[(i-1,j-1)]
        except:
            scoreMat[(i-1,j-1)] = 0

        try:
            scoreMat[(i-1,j)]
        except:
            scoreMat[(i-1,j)] = 0

        try:
            scoreMat[(i,j-1)]
        except:
            scoreMat[(i,j-1)] = 0

        try:
            scoreMat[(i,j)]
        except:
            scoreMat[(i,j)] = 0

        if scoreMat[(i,j)] == scoreMat[(i-1,j-1)] + scoring_matrix[indexChars[0]][indexChars[1]]:
            i -= 1
            j -= 1
            ind1.append(i)
            ind2.append(j)
        elif scoreMat[(i,j)] == scoreMat[(i,j-1)] + scoring_matrix[-1][indexChars[1]]:
            j -= 1
        elif scoreMat[(i,j)] == scoreMat[(i-1,j)] + scoring_matrix[indexChars[0]][-1]:
            i -= 1
        elif scoreMat[(i,j)] == 0:
            break

    ind1 = ind1[::-1]
    ind2 = ind2[::-1]

    return maxScore, ind1, ind2
