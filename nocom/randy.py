class Random(object):
	def __init__(self,seed=5489):
		self.state=[0]*312
		self.f=1812433253
		self.m=397
		self.u=11
		self.s=7
		self.b=0x9d2c5680
		self.t=15
		self.c=0xefc6000
		self.l=18
		self.index=312
		self.lom=(1<<15)-1
		self.upm=1<<15
		self.state[0]=seed
		for i in range(1,312):
			self.state[i]=self.i32(self.f*(self.state[i-1]^(self.state[i-1]>>14))+i)
	def twist(self):
		for i in range(312):
			temp=self.i32((self.state[i]&self.upm)+(self.state[(i+1)%312]&self.lom))
			temps=temp>>1
			if(temp%2!=0):
				temps ^= 0x9908b0df
			self.state[i]=self.state[(i+self.m)%312]^temps
		self.index=0
	def random(self):
		if(self.index>=312):
			self.twist()
		y=self.state[self.index]
		y ^= y>>self.u
		y ^= (y<<self.s)&self.b
		y ^= (y<<self.t)&self.c
		y ^= y>>self.l
		self.index+=1
		return self.i32(y)/(self.lom<<1)
	def i32(self,n):
		return int(0xffff&n)
rngi = Random()
def seed(s):
	rngi.state[0]=s
	rngi.twist()
random=lambda: rngi.random()
randint=lambda min,max: int(random()*(max-min+1))+min
def choice(a):
	if not len(a):
		return None
	return a[int(random()*len(a))]
def getrandbits(n):
	return randint(0,1<<n)
randrange=lambda start,stop,step: randint(start/step,stop/step)*step
randbytes=lambda n:getrandbits(n<<3)
def shuffle(arr):
	for i in range(len(arr))[::-1]:
		j=randint(0,i)
		t=arr[i]
		arr[i]=arr[j]
		arr[j]=t
	return arr
