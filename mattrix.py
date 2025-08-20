import matthy
def dim(mat):
  if(len(mat.v)==0 or len(mat.v[0]) == 0):
    return [0,0]
  return [len(mat.v),len(mat.v[0])]
mean=lambda v: sum(v)/len(v) # quick mean calculation
zeros=lambda m,n: M([[0]*n for i in range(m)]) # m x n matrix filled with zeros
ones=lambda m,n: M([[1]*n for i in range(m)]) # m x n matrix filled with ones
eye=lambda m,n: M([([0]*i)+[1]+([0]*(n-i-1)) for i in range(m)]) # m x n matrix with diagonals=1 and everything else=0
R=lambda theta: M([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]]) # 2D rotation matrix generator, theta in radians
Rx=lambda theta: M([[1,0,0],[0,cos(theta),-sin(theta)],[0,sin(theta),cos(theta)]]) # 3D x-axis rotation matrix generator
Ry=lambda theta: M([[cos(theta),0,sin(theta)],[0,1,0],[-sin(theta),0,cos(theta)]]) # 3D y-axis rotation matrix generator
Rz=lambda theta: M([[cos(theta),-sin(theta),0],[sin(theta),cos(theta),0],[0,0,1]]) # 3D z-axis rotation matrix generator
gen=lambda f,r,c: M([[f(ri,ci) for ci in range(c)] for ri in range(r)]) # Generate r x c matrix based on the current row and column
class M: # Matrix class
    def __init__(self,v): # v: number[][]
        self.m=len(v)
        self.n=len(v[0])
        self.v=v
    def apply(self,f): # apply a function to each row and column's value. Function f takes the current index value, row index, and column index
        for ri in range(len(self.v)):
            r=self.v[ri] # row vals
            for ci in range(len(r)):
                vv=r[ci] # index val
                self.v[ri][ci]=f(vv,ri,ci) # apply f
        return self
    def reduce(self,f): # f should take a 1D array 
        return f([f(r) for r in self.v])
    def kmul(self,k): # scalar multiplication of matrix
        return M([[k*v for v in il] for il in self.v])
    def sum(self): # add all the values in the matrix
        return self.reduce(sum)
    def min(self): # get the min value in the matrix
        return self.reduce(min)
    def max(self): # get the max value in the matrix
        return self.reduce(max)
    def mul(self,mat): # multiply an a x b and b x c matrix
        m0=dim(self)[0] # get matrix dimensions for self
        n0=dim(self)[1]
        m1=dim(mat)[0] # get matrix dimensions for other matrix "mat"
        n1=dim(mat)[1]
        if m1!=n0: # if dimensions aren't correct for matrix multiplication, raise error
            raise ValueError("Dimensions "+str(m0)+"x"+str(n0)+" and "+str(m1)+"x"+str(n1)+" do not match for multiplication")
        res=zeros(m0,n1) # make result matrix (a x c)
        for m in range(m0):
            for n in range(n1):
                for p in range(n0):
                    res.v[m][n]+=self.v[m][p]*mat.v[p][n]
        return res
    def add(self,mat): # add 2 m x n matrices
        if not dim(self)==dim(mat): # if dims don't match, raise error
            raise ValueError("Dimensions "+str(dim(self)[0])+"x"+str(dim(self)[1])+" and "+str(dim(mat)[0])+"x"+str(dim(mat)[1])+" do not match for addition")
        return M([[mat.v[i][j]+self.v[i][j] for j in range(len(self.v[i]))] for i in range(len(self.v))]) # do the addition in one line
    def neg(self): # negate matrix
        return kmul(self,-1)
    def sub(self,mat): # subtract 2 matrices based on .add(...)
        return self.add(mat.kmul(-1)) # just do add but multiply other matrix "mat" by a scalar -1
    def det(self): # calculate determinant
        sz=dim(self)
        if(sz[0]!=sz[1] or 0 in sz): # need square matrix nonzero-sized for determinant
            raise ValueError("Matrix has no determinant")
        bsz=dim(self)[0]
        mat=self.v # get the raw number[][] values
        if(bsz==1): # if it's a 1x1 matrix, the determinant is just the one value in it
            return mat[0][0]
        if(bsz==2): # if it's a 2x2 matriux, the values are ad-bc 
            return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]
        return sum([((-1)**c)*mat[0][c]*minor(self,0,c).det() for c in range(bsz)]) # calculate the determinant based on minors otherwise
    def tp(self): # transpose matrix (turn 90deg)
        return M(list(map(list,zip(*self.v))))
    def hdmd(self,mat): # get Hadamard product of 2 m x n matrices
        d1=dim(self)
        d2=dim(mat)
        if(d1!=d2):
            d1l=str(d1[0])+"x"+str(d1[1])
            d2l=str(d2[0])+"x"+str(d2[1])
            raise ValueError("Dimensions "+d1l+" and "+d2l+" are incompatible for Hadamard multiplication")
        return self.apply(lambda v,r,c: v*mat.v[r][c])
    def inv(self): # get the inverse of a m x m matrix
        mat=self.v
        dt=self.det()
        sz=dim(self)
        if(sz[0]!=sz[1] or dt==0): # must be a square non-zero-sized matrix
            raise ValueError("Matrix not invertible, must be a square and non-null matrix")
        szb=sz[0] # get value "m" in "m x m"
        if(szb==1):
            return M([[1/self.v[0][0]]])
        if(szb==2): # use the inversion formula for 2x2 straight up
            return M([[mat[1][1],-mat[0][1]],[-mat[1][0],mat[0][0]]]).kmul(1/dt)
        cof=[] # for bigger matrices, use the general formula for 3x3 and bigger
        for r in range(szb):
            cofr=[]
            for c in range(szb):
                mnr=_minor(mat,r,c)
                cofr.append(((-1)**(r+c))*M(mnr).det())
        cof.append(cofr)
        cof=M(cof).tp()
        return cof.kmul(1/dt)
def elim(mat): # gaussian elimination of (m) x (m+1) matrix where the last column's values are the other end of the equation 
    m = mat.v
    n=len(m)
    for i in range(n):
        if m[i][i]==0:
            return None
        for j in range(i+1,n):
            sf=m[j][i]/m[i][i]
            for k in range(n+1):
                m[j][k]-=sf*m[i][k]
    x=[0]*len(m)
    x[n-1]=m[n-1][n]/m[n-1][n-1]
    for i in range(n-2,-1,-1):
        x[i]=m[i][n]
        for j in range(i+1,n):
            x[i]=x[i]-m[i][j]*x[j]
        x[i]=x[i]/m[i][i]
    return x
def pr(pts,k): # calculate polynomial regression for k-degreed polynomial with pts in format of number[][2]
    mat=[]
    for row in range(k):
        mat.append([])
        for col in range(k):
            mat[row].append(0)
            for i in range(len(pts)):
                mat[row][col]+=pts[i][0]**(row+col)
        mat[row].append(0)
        for i in range(len(pts)):
            mat[row][-1]+=pts[i][1]*pts[i][0]**row
    return elim(M(mat))
def prfunc(pts,k): # use pr(...) but make into actual py function
    betas=pr(pts,k)
    def f(x):
        v=0
        for i in range(k):
            v+=betas[i]*x**i
        return v
    return f
def cod(pts,f): # calculate R^2 value (coefficient of determination)
    rss=sum([(pt[1]-f(pt[0]))**2 for pt in pts]) # residual sum of squares
    yb=mean([pt[1] for pt in pts]) # get y-bar (mean of y values)
    tss=sum([(pt[1]-yb)**2 for pt in pts]) # total sum of squares
    return 1-rss/tss # calculate R^2 from these vals
def _minor(mat,i,j): # raw 2D array minor 
    return [row[:j]+row[j+1:] for row in (mat[:i]+mat[i+1:])]
def minor(mat,i,j): # matrix class-based minor
    return M(_minor(mat.v,i,j))    
