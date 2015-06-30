"""
Implementation of the Kohonen Self-Similarity Map [1].

[1] Kohonen, Teuvo (1982). "Self-Organized Formation of Topologically Correct Feature Maps". Biological Cybernetics 43 (1): 59-69
"""
import numpy as np
from math import *
from sklearn.decomposition import PCA

pi = np.pi

def mag(v):
	""" returns magnitude of vector """
	return sqrt(np.dot(v,v))

def EuclideanDistance(v1, v2):
	""" usual Euclidean distance metric"""
	return mag(v1 - v2)

def MinEuclidean(v1, v2):
	""" finds the minimum (Euclidean) distance between two vectors, 
		allowing for an offset between their indices. (think phase shift...)
	"""
	mind = None
	for i in range(len(v1)):
		v1p = [ v1[j - i] for j in range(len(v1)) ]
		d = EuclideanDistance(v1p, v2)
		if mind is None or d < mind:
			mind = d
	return mind

gauss = lambda x : exp(- x * x/2.)

def GaussianTheta(N, Nsteps):
	""" Returns a Theta function that is one at the BMU index, and
		decreases as ~e(-d^2/2s^2) away from the BMU index; 
		sigma is initially equal to the number of grid points along one dimension,
	"""
	def func(u, v, s):
		D = np.array(u) - np.array(v)
		#print D
		U = sqrt(np.dot(D, D))/N
		dx = 1./N
		sigma = 1 - (1 - dx)*(float(s)/float(Nsteps))
		return gauss(U/sigma)
	return func


class KohonenMap:
	def __init__(self, metric_function=EuclideanDistance,theta=None, alpha=None, N=25, Nsteps=10):
		#if Nd not in [ 1, 2, 3 ] : raise ValueError("Nd must be 1, 2, or 3.")
		#self.shape = tuple([ N for i in range(Nd) ]) 
		self.N = N
		self.shape = (N, N)
		self.Nsteps = Nsteps
		
		if not theta is None:
			self.Theta = theta
		else:
			self.Theta = GaussianTheta(N, Nsteps)

		if not alpha is None:
			self.Alpha = alpha
		else:
			# linear weight function...
			self.Alpha = lambda s : 1. - 0.9*float(s)/float(Nsteps)

	def _init_grid(self, observations, initialization_type="PCA"):


		self.grid = []
		data_size = len(observations[0])
		if initialization_type == "random":
			for i in range(self.N):
				self.grid.append([ ])
				for j in range(self.N):
					V = np.random.random(data_size)
					self.grid[i].append(V/mag(V))
		elif initialization_type == "PCA":
			P = PCA(n_components=2).fit(observations).components_
			v1 = P[0]/mag(P[0])
			v2 = P[1]/mag(P[1])
			s = np.linspace(0,1,self.N)
			for i,a in enumerate(s):
				self.grid.append([ ])
				for j,b in enumerate(s):
					v = a*v1 + b*v2
					if mag(v) > 0: v /= mag(v)
					self.grid[i].append(v)

		else:
			raise ValueError("Do not recognize initialization_type")

	def _find_bmu(self, D):
		""" 
		Slow, naive search for closest match in the K map. 
		"""
		mind, bmui = None, None
		for i in range(self.N):
			for j in range(self.N):
				vd = D - self.grid[i][j]
				d = mag(vd)
				
				if (mind is None) or (d < mind):
					mind = d
					bmui = (i, j)
		return bmui

	def fit(self, observations,**kwargs):

		self._init_grid(observations, **kwargs)

		# Naive training. Maybe we can speed this up by
		# 	* doing a more heirarchical search (Ng^2 -> NglogNg)
		#		* user specifies minN, which indicates the top-most
		#		  level of the training.
		#	* stochastically estimating (D - Wvs)
		for s in range(self.Nsteps):
			print "Training step ",s," of ", self.Nsteps
			for D in observations:
				
				# find BMU
				bmui = self._find_bmu(D)
				#print bmui

				# Imprint BMU on grid.
				for i in range(self.N):
					for j in range(self.N):
						Wvs = self.grid[i][j]
						Wvss = Wvs + self.Theta(bmui, (i,j), s) * self.Alpha(s) * (D - Wvs)
						self.grid[i][j] = Wvss


	def transform(self, observations, **kwargs):
		results = [ ]
		for D in observations:
			u, v = self._find_bmu(D)
			results.append(np.array([ float(u)/self.N, float(v)/self.N ]))
		return results
