from matthy import *
from matthext import *
from mattrix import *
from mattvec import *
from randy import *
def mat2vec(mat):
		if(len(mat.v)==1):
				return V(mat.v[0])
		return [V(arr) for arr in mat.v]
mat2vec_ttb = lambda mat: mat2vec(mat.tp())
mat2vec_ltr = mat2vec
