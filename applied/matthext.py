import matthy
from matthy import sigma
addsign=lambda x: ("+" if x>=0 else "")+str(x)
def lim(appr,f,d=7):
  xvl=[appr-1]
  yvl=[f(xvl[0])]
  xvr=[appr+10**-d]
  yvr=[f(xvr[0])]
  ml=[]
  mr=[]
  for i in range(1,d):
    hl=10**(-i)
    hr=10**(i-d)
    xvl.append(appr-hl)
    yvl.append(f(appr-hl))
    xvr.append(appr+hr)
    yvr.append(f(appr+hr))
    ml.append((yvl[i]-yvl[i-1])/(xvl[i]-xvl[i-1]))
    mr.append((yvr[i]-yvr[i-1])/(xvr[i]-xvr[i-1]))
  mlt=trend(ml)
  mrt=trend(mr)
  return [yvr[0],mlt,mrt]
def trend(m,d=0.1):
  pv=m[0]
  ups=0
  downs=0
  for v in m[1:]:	 
    if(v-pv>d):
      ups+=1
    if(v-pv<-d):
      downs+=1
    pv=v
  return (ups-downs)/len(m)
def pcert(ret):
  lav=ret[1]
  rav=ret[2]
  tab = ["incr./decr.","incr.","decr."]
  f = lambda v: tab[matthy.sign(v)]
  print(str(abs(int(lav*100)))+"% of delta-slopes from left were "+f(lav))
  print(str(abs(int(rav*100)))+"% of delta-slopes from right were "+f(rav))
  print("a close value calculated to "+str(ret[0]))
def ddx(x,f):
  return lim(0,lambda h: (f(x+h)-f(x))/h)
def slr(pts):
  dx=sum([pt[0] for pt in pts])/len(pts)
  dy=sum([pt[1] for pt in pts])/len(pts)
  bh=sigma(0,len(pts)-1,lambda i: (pts[i][0]-dx)*(pts[i][1]-dy))
  bh/=sigma(0,len(pts)-1,lambda i: (pts[i][0]-dx)**2)
  ah=dy-bh*dx
  return [bh,ah,"y="+str(bh)+"x"+addsign(ah)]
