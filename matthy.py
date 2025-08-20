accuracy=25
liaccuracy=accuracy*100
infim=32767
pi=3.141592653589793
tau=2*pi
e=2.718281828459045
ggam=5
ngam=5
pgam=[
	1.0000018972739440364,
	76.180082222642137322,
	-86.505092037054859197,
	24.0128985819226859,
	-1.2296028490285820771
]
def factorial(n):
	if(n==0): 
		return 1
	r = 1
	for i in range(1,n+1):
		r*=i
	return r
def comb(n,k):
	return factorial(n)/(factorial(k)*factorial(n-k))
def trunc(x):
	return x-(x%1)
def floor(x):
	if(x==trunc(x)):
		return x
	return trunc(x)+1
trunc=floor
def sign(x):
	if(x==0):
		return 0
	elif(x>0):
		return 1
	else:
		return -1
def copysign(x,y):
	return x*sign(y)
# sigma(...) is NOT in default math
def sigma(min,max,f):
	r = 0
	for i in range(min,max+1):
		r += f(i)
	return r
def exp(z):
	return sigma(0,accuracy,lambda n: z**n/factorial(n))
def ln(v,acc=accuracy):
	x = v
	y = 1
	for i in range(1,accuracy):
		y += 2*((x-exp(y))/(x+exp(y)))
	return y
def log(x,b=exp(1)):
	return ln(x)/ln(b)
log1p = lambda x: ln(1+x)
log2 = lambda x: log(x,2)
log10 = lambda x: log(x,10)
pow = lambda x,y: x**y
exp2 = lambda x: 2**x
expm1 = lambda x: exp(x)-1
sqrt = lambda x: x**0.5
cbrt = lambda x: x**(1/3)
degrees = lambda x: x*180/pi
radians = lambda x: x*pi/180
cos = lambda x: sigma(0,accuracy,lambda n: (-1)**n*x**(2*n)/factorial(2*n))
sin = lambda x: sigma(0,accuracy,lambda n: (-1)**n*x**(2*n+1)/factorial(2*n+1))
tan = lambda x: sin(x)/cos(x)
sec = lambda x: 1/cos(x)
csc = lambda x: 1/sin(x)
cot = lambda x: 1/tan(x)
cosh = lambda x: (exp(x)+exp(-x))/2
sinh = lambda x: (exp(x)-exp(-x))/2
tanh = lambda x: sinh(x)/cosh(x)
sech = lambda x: 1/cosh(x)
csch = lambda x: 1/sinh(x)
coth = lambda x: 1/tanh(x)
asin = lambda x: sigma(0,accuracy,lambda k: asins(k,x))
asins = lambda k,x: comb(2*k,k)/((4**k)*(2*k+1))*x**(2*k+1)
acos = lambda x: pi/2-asin(x)
atan = lambda x,acc=accuracy: sigma(0,acc,lambda n: (-1)**n*(x**(2*n+1)/(2*n+1)))
atan2 = lambda y,x,acc=accuracy: atan(y/x,acc)
asinh = lambda x: ln(x+sqrt(x**2+1))
acosh = lambda x: ln(x+sqrt(x**2-1))
atanh = lambda x: ln((1+x)/(1-x))/2
acsch = lambda x: ln(1/x+sqrt(1/x**2+1))
asech = lambda x: ln(1/x+sqrt(1/x**2-1))
acoth = lambda x: ln((x+1)/(x-1))
erf = lambda x: 2/sqrt(pi)*sigma(0,accuracy,lambda n: ((-1)**n*x**(2*n+1))/((2*n+1)*factorial(n)))
erfc = lambda x: 1-erf(x)
def perm(n,k=None):
	if(k==None):
		k=n
	if(k<=n):
		return factorial(n)/factorial(n-k)
	return 0
isclose=lambda a,b,rel_tol=1e-9: abs(a-b)<rel_tol
fabs=lambda x: sign(x)*x
fma=lambda x,y,z: (x*y)+z
fmod=lambda x,y:x-trunc(x/y)*y
ldexp=lambda x,i:x*2**i
dist=lambda p,q: sqrt(sum((px-qx)**2 for px,qx in zip(p,q)))
hypot=lambda coordinates: sqrt(sum([x**2 for x in coordinates]))

def riemann(a,b,n,f):
	return sigma(1,n,lambda i: ((b-a)/n)*f(a+i*((b-a)/n)))
def capi(min,max,f,sv=1):
	s=sv
	for i in range(min,max+1):
		s*=f(i)
	return s
def gamma(z,acc=liaccuracy):
	return capi(1,acc,lambda n:(1+1/n)**z/(1+z/n),sv=1/z)
def gamma2(z):
	p=pgam
	n=ngam
	g=ggam
	if(z<0.5):
		return pi/(sin(pi*z)*gamma(1-z))
	z-=1
	x=p[0]
	x+=sigma(1,4,lambda i: p[i]/(z+i))
	t=z+g+0.5
	return sqrt(tau)*t**(z+0.5)*exp(-t)*x	 
	return riemann(0,32767,50,lambda t: itgamma(t,z))
def agm(x,y):
	if(not (x>=y and y>=0)):
		raise ValueError("x and y must fit x>=y>=0")
	a=x
	b=y
	for i in range(accuracy):
		at=a
		a=(at+b)/2
		b=sqrt(at*b)
	return a
gm = lambda x,y: sqrt(x*y) 
G = 1/agm(sqrt(2),1)
lem = pi*G
def hm(*a):
	return len(a)/sigma(0,len(a)-1,lambda x: 1/a[x])
zeta = lambda s: sigma(1,liaccuracy,lambda n: n**-s)
