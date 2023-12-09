# 1. Metode eliminasi Gauss-Jordan
def MakeIdentity(n):
    # KAMUS LOKAL
    # M: array of array of float
    # i, j: integer
    # ALGORITMA 
    M = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if (i == j):
                M[i][j] = 1
            else:
                M[i][j] = 0
    return M

def SwapRow(M,a,b):
    # KAMUS LOKAL
    # temp: array of float
    # ALGORITMA
    temp = [0 for i in range(len(M[0]))]
    temp = M[a]
    M[a] = M[b]
    M[b] = temp

def MultiplyRow(M,row,c):
    # KAMUS LOKAL
    # i: integer
    # ALGORITMA
    for i in range(len(M[0])):
        M[row][i] *= c

def AddRow(M,row1,row2,c):
    # KAMUS LOKAL
    # i: integer
    # ALGORITMA
    for i in range(len(M)):
        M[row2][i] += (M[row1][i] * c)

def InversGaussJordan(Min):
    # KAMUS LOKAL
    # Mout: array of array of float
    # i, j: integer
    # ALGORITMA
    Mout = MakeIdentity(len(Min))
    for i in range(len(Min)):
        if (Min[i][i] == 0):
            j = i+1
            while (Min[j][i] == 0):
                j += 1
            SwapRow(Min,i,j)
            SwapRow(Mout,i,j)

        MultiplyRow(Mout,i,1/Min[i][i])
        MultiplyRow(Min,i,1/Min[i][i])

        for j in range(len(Min)):
            if (j != i):
                AddRow(Mout,i,j,(-1)*Min[j][i])
                AddRow(Min,i,j,(-1)*Min[j][i])
    return Mout

# 2. Metode matriks adjoin
def Det(M):
    # KAMUS LOKAL
    # Ctr, i, j: integer
    # D: float
    # ALGORITMA
    Ctr = 1
    for i in range(len(M)):
        if (M[i][i] == 0):
            Ctr *= (-1)
            j = i + 1
            while (M[j][i] == 0):
                j += 1
            SwapRow(M,i,j)
        for j in range(len(M)):
            AddRow(M,i,j,(-1)*(M[j][i]/M[i][i]))
    D = 1
    for i in range(len(M)):
        D *= M[i][i]
    return (Ctr*D)

def MinorMatrix(M,row,col):
    # KAMUS LOKAL
    # Mo: array of array of float
    # i, j, k, l: integer
    # ALGORITMA
    Mo = [[0 for i in range(len(M)-1)] for i in range(len(M)-1)]
    l = 0
    for i in range (len(M)):
        if (i != row):
            k = 0
            for j in range(len(M)):
                if (j != col):
                    Mo[l][k] = M[i][j]
                    k += 1
            l += 1
    return Mo

def CofactorMatrix(M):
    # KAMUS LOKAL
    # C: array of array of float
    # i, j: integer
    # Ctr, Temp: integer
    # ALGORITMA
    C = [[0 for i in range(len(M))] for i in range(len(M))]
    Ctr = -1
    for i in range(len(M)):
        Ctr *= -1
        Temp = Ctr
        for j in range(len(M)):
            C[i][j] = Temp*Det(MinorMatrix(M,i,j))
            Temp *= -1
    return C

def Transpose(M):
    # KAMUS LOKAL
    # i, j: integer
    # Tm: array of array of float
    # ALGORITMA
    Tm = [[0 for i in range(len(M))] for i in range(len(M))]
    for i in range(len(M)):
        for j in range(len(M)):
            Tm[i][j] = M[j][i]
    return Tm

def InversAdjoin(Min):
    # KAMUS LOKAL
    # Mout: array of array of float
    # D: float
    # i: integer
    # ALGORITMA
    Mout = CofactorMatrix(Min)
    Mout = Transpose(Mout)
    D = Det(Min)
    for i in range(len(Min)):
        MultiplyRow(Mout,i,1/D)
    return Mout