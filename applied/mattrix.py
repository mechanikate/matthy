import matthy
def dim(mat):
  if(len(mat.v)==0 or len(mat.v[0]) == 0):
    return [0,0]
  return [len(mat.v),len(mat.v[0])]
def dimt(mat):
    d = dim(mat)
    return (d[0],d[1])
mean=lambda v: sum(v)/len(v) 
zeros=lambda m,n: M([[0]*n for i in range(m)]) 
ones=lambda m,n: M([[1]*n for i in range(m)]) 
eye=lambda m,n: M([([0]*i)+[1]+([0]*(n-i-1)) for i in range(m)]) 
R=lambda theta: M([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]]) 
Rx=lambda theta: M([[1,0,0],[0,cos(theta),-sin(theta)],[0,sin(theta),cos(theta)]]) 
Ry=lambda theta: M([[cos(theta),0,sin(theta)],[0,1,0],[-sin(theta),0,cos(theta)]]) 
Rz=lambda theta: M([[cos(theta),-sin(theta),0],[sin(theta),cos(theta),0],[0,0,1]]) 
gen=lambda f,r,c: M([[f(ri,ci) for ci in range(c)] for ri in range(r)]) 
class M: 
    def __init__(self,v): 
        self.m=len(v)
        self.n=len(v[0])
        self.v=v
    def apply(self,f): 
        for ri in range(len(self.v)):
            r=self.v[ri] 
            for ci in range(len(r)):
                vv=r[ci] 
                self.v[ri][ci]=f(vv,ri,ci) 
        return self
    def reduce(self,f): 
        return f([f(r) for r in self.v])
    def kmul(self,k): 
        return M([[k*v for v in il] for il in self.v])
    def sum(self): 
        return self.reduce(sum)
    def min(self): 
        return self.reduce(min)
    def max(self): 
        return self.reduce(max)
    def mul(self,mat): 
        m0=dim(self)[0] 
        n0=dim(self)[1]
        m1=dim(mat)[0] 
        n1=dim(mat)[1]
        if m1!=n0: 
            raise ValueError("Dimensions "+str(m0)+"x"+str(n0)+" and "+str(m1)+"x"+str(n1)+" do not match for multiplication")
        res=zeros(m0,n1) 
        for m in range(m0):
            for n in range(n1):
                for p in range(n0):
                    res.v[m][n]+=self.v[m][p]*mat.v[p][n]
        return res
    def add(self,mat): 
        if not dim(self)==dim(mat): 
            raise ValueError("Dimensions "+str(dim(self)[0])+"x"+str(dim(self)[1])+" and "+str(dim(mat)[0])+"x"+str(dim(mat)[1])+" do not match for addition")
        return M([[mat.v[i][j]+self.v[i][j] for j in range(len(self.v[i]))] for i in range(len(self.v))]) 
    def neg(self): 
        return kmul(self,-1)
    def sub(self,mat): 
        return self.add(mat.kmul(-1)) 
    def det(self): 
        sz=dim(self)
        if(sz[0]!=sz[1] or 0 in sz): 
            raise ValueError("Matrix has no determinant")
        bsz=dim(self)[0]
        mat=self.v 
        if(bsz==1): 
            return mat[0][0]
        if(bsz==2): 
            return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]
        return sum([((-1)**c)*mat[0][c]*minor(self,0,c).det() for c in range(bsz)]) 
    def tp(self): 
        self = M(list(map(list,zip(*self.v))))
        return self
    def hdmd(self,mat): 
        d1=dim(self)
        d2=dim(mat)
        if(d1!=d2):
            d1l=str(d1[0])+"x"+str(d1[1])
            d2l=str(d2[0])+"x"+str(d2[1])
            raise ValueError("Dimensions "+d1l+" and "+d2l+" are incompatible for Hadamard multiplication")
        return self.apply(lambda v,r,c: v*mat.v[r][c])
    def inv(self): 
        mat=self.v
        dt=self.det()
        sz=dim(self)
        if(sz[0]!=sz[1] or dt==0): 
            raise ValueError("Matrix not invertible, must be a square and non-null matrix")
        szb=sz[0] 
        if(szb==1):
            return M([[1/self.v[0][0]]])
        if(szb==2): 
            return M([[mat[1][1],-mat[0][1]],[-mat[1][0],mat[0][0]]]).kmul(1/dt)
        cof=[] 
        for r in range(szb):
            cofr=[]
            for c in range(szb):
                mnr=_minor(mat,r,c)
                cofr.append(((-1)**(r+c))*M(mnr).det())
        cof.append(cofr)
        cof=M(cof).tp()
        return cof.kmul(1/dt)
    def off_diag_frobnorm(self):
        m,n=dimt(self)
        if(m!=n):
            raise ValueError("Matrix must be square")
        return sum([sum([self.v[i][j]**2 if i != j else 0 for j in range(n)]) for i in range(n)])**0.5
    def eigvals(self,accuracy=matthy.accuracy,tol=1e-9,max_iters=1e4):
        mat = M(self.v)
        sized_eye=eye(*dimt(self))
        frob_norm = matthy.infim
        iters = 0
        while frob_norm > tol and iters < max_iters:
            mu=sized_eye.kmul(mat.v[-1][-1]+tol/10)
            Q,R=gramschmidt(mat.sub(mu))
            mat=R.mul(Q).add(mu)
            frob_norm=mat.off_diag_frobnorm()
            iters += 1
        return [mat.v[i][i] for i in range(min(dim(mat)))]

def elim(mat): 
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
def pr(pts,k): 
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
def prfunc(pts,k): 
    betas=pr(pts,k)
    def f(x):
        v=0
        for i in range(k):
            v+=betas[i]*x**i
        return v
    return f
def cod(pts,f): 
    rss=sum([(pt[1]-f(pt[0]))**2 for pt in pts]) 
    yb=mean([pt[1] for pt in pts]) 
    tss=sum([(pt[1]-yb)**2 for pt in pts]) 
    return 1-rss/tss 
def _minor(mat,i,j): 
    return [row[:j]+row[j+1:] for row in (mat[:i]+mat[i+1:])]
def minor(mat,i,j): 
    return M(_minor(mat.v,i,j))   
def gramschmidt(mat):
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
    return (Q,R)

