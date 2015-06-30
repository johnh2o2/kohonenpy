import numpy as np
from math import *
from kohonen import *
import matplotlib.pyplot as plt

Nobs = 1000
Nh = 3
Nt = 4
Nbs = 20
Ngrid = 15
Nsteps = 3*Ngrid
types = []
for t in range(Nt):
	types.append({ 'a' : [  np.random.random() for i in range(Nh) ], 
					'ph' : [  2*pi*np.random.random() for i in range(Nh) ] })

def make_lc(info):
	w = np.random.random()
	ph = 2*pi*np.random.random()
	sigma = 0.05
	phs = np.linspace(0,1)

	fluxes = np.array([ sum([ 0.5 * (A*cos((i + 1)*2*pi*PHASE - ph - PH) + 1) for i, A, PH in zip(np.arange(len(info['a'])),info['a'], info['ph']) ]) + 0.05*np.random.normal() for PHASE in phs])


	binned_flux = np.zeros(Nbs)
	binned_phases = np.linspace(0,1, Nbs)
	binned_phases += binned_phases[0]
	numbers = np.zeros(Nbs)
	for PH, FLUX in zip(phs, fluxes):
		i=0
		while PH > binned_phases[i]: i+=1
		binned_flux[i] += FLUX
		numbers[i] += 1.
	binned_flux/=numbers
	return binned_phases, binned_flux

Map = KohonenMap(N=Ngrid, Nsteps=Nsteps, metric_function=MinEuclidean)
Nobs_per_type = Nobs/Nt
observations = []
all_obs = []
for Type in types:
	type_obs = []
	for i in range(Nobs_per_type):
		binned_phases, binned_fluxes = make_lc(Type)
		Fnorm = mag(binned_fluxes)

		observations.append(binned_fluxes/Fnorm)
		type_obs.append(binned_fluxes/Fnorm)
	all_obs.append(type_obs)
np.random.shuffle(observations)
Map.fit(observations, initialization_type="random")

#P = PCA(n_components=2).fit(observations).components_
#ax.plot(binned_phases, P[0])
#ax.plot(binned_phases, P[1])
#ax.plot(binned_phases, Map.grid[0][0])
#ax.plot(binned_phases, Map.grid[1][1])
#ax.plot(binned_phases, Map.grid[3][5])
#ax.plot(binned_phases, P[2])
#plt.show()

colors = [ 'r', 'b', 'c', 'g' ]


for i, Type in enumerate(types):
	f = plt.figure()
	ax = f.add_subplot(111)
	#print all_obs[i]
	results = Map.transform(all_obs[i])
	grid = np.zeros((Ngrid,Ngrid))
	x = np.linspace(0, 1, Ngrid)
	#print results 
	v1s = [ r[0] for r in results ]
	v2s = [ r[1] for r in results ]
	dx = x[1] - x[0]
	for X, Y in zip(v1s, v2s):
		grid[int(X/dx)][int(Y/dx)]+= 1./len(all_obs[i])

	c = colors[i]
	#ax.scatter(v1s, v2s, marker='o', facecolor=c, alpha=0.2,edgecolor=None,linewidth=0.5)
	heatmap = ax.pcolor(grid, cmap=plt.cm.Blues)
	plt.colorbar(heatmap)
	ax.set_title("type %d"%(i))
plt.show()









