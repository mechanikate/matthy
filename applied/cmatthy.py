import matthy
from matthy import pi, isclose
def shiftdomain(v,mi1,ma1,mi2,ma2):
  if(v==mi1):
    return mi2
  if(v==ma1):
    return ma2
  return (((v-mi1)*(ma2-mi2))/(ma1-mi1))+mi2
testv = [
  [0,0,0],
  [1,1,pi/4],
  [0,1,pi/2],
  [-1,1,3*pi/4],
  [-1,0,pi],
  [-1,1,3*pi/4],
  [0,-1,-pi/2],
  [1,-1,-pi/4]
]
def ctest(f,v,rel=1e-4):
  passed=0
  failed=0
  total=len(v)
  for i in range(len(v)):
    e=v[i]
    res=f(complex(e[0],e[1]))
    if(isclose(res,e[2],rel)):
      passed+=1
    else:
      failed+=1
      print("fail @ "+str(i)+": f("+str(e[0])+"+"+str(e[1])+"i = "+str(e[2])+" != "+str(res))
  print(str(passed)+"/"+str(total))
  
class complex:
  def __init__(self,real,imag):
    self.real = real
    self.imag = imag
  def abs(self):
    return matthy.sqrt(self.real**2+self.imag**2)

def csqrt(z):
  az = z.abs()
  den = complex(az+z.real,z.imag).abs()
  re = z.real+az
  im = z.imag
  re*=matthy.sqrt(az)/den
  im*=matthy.sqrt(az)/den
  return complex(re,im)

def phase(x,acc=1500):
  if(x.real==0):
    if(x.imag==0):
      return 0
    return pi/2 if x.imag >= 0 else -pi/2
  r = matthy.atan2(x.imag,x.real,acc)
  return r
  
polar = lambda x: (x.abs(), phase(x))
rect = lambda r,phi: complex(r*matthy.cos(phi),r*matthy.sin(phi))
