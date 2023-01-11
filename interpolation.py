from numpy import linspace, array

# Interpolation bilinéaire

def lineaire(L,gros):
    T = []
    for i in range(len(L)-1):
        for t in linspace(0,1,gros + 1)[:-1]:
            T.append((1-t)*L[i] + t*L[i+1])
    return T + [L[-1]]

def transpose(M):
    T=[]
    for j in range(len(M[1])):
        C=[]
        for i in range(len(M)):
            C.append(M[i][j])
        T.append(C)
    return T

def bilineaire(M,gros):
    M = transpose(M)
    for i in range(len(M)):
        M[i] = lineaire(M[i],gros)
    T = transpose(M)
    for i in range(len(T)):
        T[i] = lineaire(T[i],gros)
    return T

# Interpolation par inverse des distances

def inverse(M,gros):
    n = len(M)
    T = []
    for ti in range(n*gros):
        Ti = []
        for tj in range(n*gros):
            Tij = 0
            ht = 0
            for mi in range(n):
                for mj in range(n):
                    dij = (ti - gros*mi - gros/2)**2 + (tj - gros*mj - gros/2)**2 + 20
                    Tij += M[mi][mj] / dij
                    ht += 1 / dij
            Ti.append(Tij / ht)
        T.append(Ti)
    return T

# Interpolation par voisins naturels

def sarrus(x1,x2,x3,x4,x5,x6,x7,x8,x9):
    return x1*x5*x9 + x4*x8*x3 + x7*x2*x6 - x3*x5*x7 - x6*x8*x1 - x9*x2*x4

def coordonnee(A,B,C,coord):
    xA, yA = A
    xB, yB = B
    xC, yC = C
    delta = 2 * sarrus(xA,yA,1,xB,yB,1,xC,yC,1)
    if coord == "x":
        s = sarrus(xA**2+yA**2,yA,1,xB**2+yB**2,yB,1,xC**2+yC**2,yC,1)
        return s / delta
    else:
        s = sarrus(xA**2+yA**2,xA,1,xB**2+yB**2,xB,1,xC**2+yC**2,xC,1)
        return -1 * (s / delta)

def permutation(ordre,L):
    M = []
    for i in ordre:
        M.append(L[i])
    return array(M)

def hypothenuse(ab,ac):
    return (ab**2 + ac**2)**.5

def pondSibson(gros):
    G = 2*gros
    P = []
    for _ in range(G):
        P.append([[]]*G)
    for i in range(gros):
        for j in range(i,gros):
            M = (i + .5, j + .5)
            c0 = coordonnee((0,0),(G,0),M,"y")
            c1 = coordonnee((G,0),(G,G),M,"x")
            c2 = coordonnee((0,G),(G,G),M,"y")
            c3 = coordonnee((0,0),(0,G),M,"x")
            d0 = gros - c0
            d1 = c1 - gros
            d2 = c2 - gros
            d3 = gros - c3
            A = array([d3*d0,d0*d1,d1*d2,d2*d3])
            A = A / sum(A)
            P[i][j] = permutation([0,1,2,3],A)
            P[j][i] = permutation([0,3,2,1],A)
            P[G-1-i][j] = permutation([1,0,3,2],A)
            P[G-1-j][i] = permutation([3,0,1,2],A)
            P[j][G-1-i] = permutation([1,2,3,0],A)
            P[i][G-1-j] = permutation([3,2,1,0],A)
            P[G-1-i][G-1-j] = permutation([2,3,0,1],A)
            P[G-1-j][G-1-i] = permutation([2,1,0,3],A)
    return P

def pondLaplace(gros):
    G = 2*gros
    P = []
    for _ in range(G):
        P.append([[]]*G)
    for i in range(gros):
        for j in range(i,gros):
            M = (i + .5, j + .5)
            c0 = coordonnee((0,0),(G,0),M,"y")
            c1 = coordonnee((G,0),(G,G),M,"x")
            c2 = coordonnee((0,G),(G,G),M,"y")
            c3 = coordonnee((0,0),(0,G),M,"x")
            i0 = hypothenuse(gros - c3, gros - c0)
            i1 = hypothenuse(gros - c0, c1 - gros)
            i2 = hypothenuse(c1 - gros, c2 - gros)
            i3 = hypothenuse(c2 - gros, gros - c3)
            d0 = hypothenuse(M[0],M[1])
            d1 = hypothenuse(G-M[0],M[1])
            d2 = hypothenuse(G-M[0],G-M[1])
            d3 = hypothenuse(M[0],G-M[1])
            A = array([i0/d0,i1/d1,i2/d2,i3/d3])
            A = A / sum(A)
            P[i][j] = permutation([0,1,2,3],A)
            P[j][i] = permutation([0,3,2,1],A)
            P[G-1-i][j] = permutation([1,0,3,2],A)
            P[G-1-j][i] = permutation([3,0,1,2],A)
            P[j][G-1-i] = permutation([1,2,3,0],A)
            P[i][G-1-j] = permutation([3,2,1,0],A)
            P[G-1-i][G-1-j] = permutation([2,3,0,1],A)
            P[G-1-j][G-1-i] = permutation([2,1,0,3],A)
    return P

def pondinverse(gros):
    G = 2*gros
    P = []
    for _ in range(G):
        P.append([[]]*G)
    for i in range(gros):
        for j in range(i,gros):
            M = (i + .5, j + .5)
            d0 = hypothenuse(M[0],M[1])
            d1 = hypothenuse(G-M[0],M[1])
            d2 = hypothenuse(G-M[0],G-M[1])
            d3 = hypothenuse(M[0],G-M[1])
            A = array([1/d0,1/d1,1/d2,1/d3])
            A = A / sum(A)
            P[i][j] = permutation([0,1,2,3],A)
            P[j][i] = permutation([0,3,2,1],A)
            P[G-1-i][j] = permutation([1,0,3,2],A)
            P[G-1-j][i] = permutation([3,0,1,2],A)
            P[j][G-1-i] = permutation([1,2,3,0],A)
            P[i][G-1-j] = permutation([3,2,1,0],A)
            P[G-1-i][G-1-j] = permutation([2,3,0,1],A)
            P[G-1-j][G-1-i] = permutation([2,1,0,3],A)
    return P

def voisin(M,P,gros):
    G = 2*gros
    T = array([[0.]*G*7]*G*7)
    for mi in range(7):
        for mj in range(7):
            for i in range(G):
                for j in range(G):
                    Tij  = M[mi][mj]*P[i][j][0]
                    Tij += M[mi+1][mj]*P[i][j][1]
                    Tij += M[mi+1][mj+1]*P[i][j][2]
                    Tij += M[mi][mj+1]*P[i][j][3]
                    T[G*mi + i][G*mj + j] = Tij
    return T
