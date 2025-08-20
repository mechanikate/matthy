from matthy import *
from matthext import *
from mattrix import *
from mattvec import *
from randy import *
def mat2vec(mat):
		if(len(mat.v)==1):
				return V(mat.v[0])
		return [V(arr) for arr in mat.v]
def dimt(mat):
    arr = dim(mat)
    return arr[0], arr[1]
def vec2mat(vecs):
        return M([vec.v for vec in vecs])
mat2vec_ttb = lambda mat: mat2vec(mat.tp())
mat2vec_ltr = mat2vec
def getcol(mat, c):
    return V([mat.v[i][c] for i in range(len(mat.v))]) 
"""def gramschmidt(mat):
    m,n=dimt(mat)
    Q=zeros(m,n)
    R=zeros(n,n)
    for j in range(n):
        Aj=[mat.v[i][j] for i in range(m)]
        for i in range(j):
            Qi = [Q.v[k][i] for k in range(m)]
            R.v[i][j] = sum(Qi[k]*Aj[k] for k in range(m))
            Aj = [Aj[k]-R.v[i][j]*Qi[k] for k in range(m)]
        R.v[j][j] = sum(x**2 for x in Aj)**0.5
        if R.v[j][j] == 0:
            raise ValueError("Columns are linearly dependent")
        for i in range(m):
            Q.v[i][j] = Aj[i]/R.v[j][j]
    return (Q,R)"""
