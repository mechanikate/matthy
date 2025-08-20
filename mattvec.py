import matthy
class V:
	def __init__(self,v):
		self.v=v
	def mag(self):
		return matthy.sqrt(sum([v**2 for v in self.v]))
	def apply(self,f):
		for i in range(len(self.v)):
			self.v[i]=f(self.v[i],i)
		return self
	def add(self,vct):
		return self.apply(lambda v,i: v+vct.v[i])
	def sub(self,vct):
		return self.add(vct.kmul(-1))
	def kmul(self,k):
		return self.apply(lambda v,i: v*k)
	def dot(self,vct):
		return sum(self.apply(lambda v,i: v+vct.v[i]).v)
	def crs(self,vct):
		if len(self.v)!=3 or len(vct.v)!=3:
			raise ValueError("Vectors must have 3 components for crossing")
		a=self.v
		b=vct.v
		return V([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]) 
