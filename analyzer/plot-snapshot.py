from math import *
import numpy as np
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colorbar as mc
#import cv2 as cv
import sys

m = 3

def Plot_File(filename1):
	global m
	f = open(filename1)

	array = [];
	flag = True

	T = 0.0

	for line in f: # read rest of lines
		if flag:
			info = line.split()
			T = float(info[2])
			frame = int(info[3])
			t = float(info[4])
			polarization = float(info[6])
			print T, frame
			flag = False
		else:
			v = [float(x) for x in line.split()]
			if (v != []):
				array.append(v)

	N = int(sqrt(len(array)))
	n = N/m

	x = np.zeros((N, N))
	y = np.zeros((N, N))
	xs = np.zeros((n, n))
	ys = np.zeros((n, n))
	density = np.zeros((N, N))
	densitys = np.zeros((n, n))
	vx = np.zeros((n, n))
	vy = np.zeros((n, n))
#	wx = np.zeros((N/n, N/n))
#	wy = np.zeros((N/n, N/n))
	for i in xrange(N):
		for j in xrange(N):
			x[N-j-1][i] = array[i*N+j][0]
			y[N-j-1][i] = array[i*N+j][1]
			density[N-j-1][i] = array[i*N+j][4]
	for i in xrange(n):
		for j in xrange(n):
			densitys[n-j-1][i] = array[m*i*N+m*j][4]
			xs[n-j-1][i] = array[m*i*N+m*j][0]
			ys[n-j-1][i] = array[m*i*N+m*j][1]
			vx[n-j-1][i] = array[m*i*N+m*j][2]
			vy[n-j-1][i] = array[m*i*N+m*j][3]


	maxc = np.amax(density)
	minc = np.amin(density)
	maxc = ceil(maxc)
	minc = floor(minc)
	midc = (maxc + minc) / 2

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif', size=20)
	#plt.rc('xtick', labelsize=20) 
	#plt.rc('ytick', labelsize=20) 


	fig = plt.figure()

	ax = fig.add_subplot()

	maxc = np.amax(density)
	minc = np.amin(density)
	maxc = ceil(maxc)
	minc = floor(minc)
	maxc = ceil(maxc)
	minc = floor(minc)
	midc = (maxc + minc) / 2
	if (maxc == minc):
		maxc = midc + 0.1
		minc = midc - 0.1

	def formater():
		if (maxc >= 100):
			return '%1.0f'
		else:
			return '%1.1f'
	plt.imshow(density, extent = (-0.5, 0.5, -0.5, 0.5), cmap="hot")
	Q = plt.quiver(xs,ys, vx, vy, densitys, linewidths=0.6, headwidth=3, headlength=2, headaxislength=2, width=0.007, edgecolors="#4444FF", scale=30, cmap='hot', clim=[minc,maxc])

	cb = plt.colorbar(orientation='vertical', fraction=0.082, shrink=1.0, aspect=10, pad=0.01, ticks=[minc,midc,maxc], format=formater())
	plt.clim(minc,maxc)
#plt.colorbar.make_axes(parents, location="top")
#	cb.set_label(r"Density", labelpad=-62)
	mytext = 'Density'+ r'$(\rho)$'
	plt.text(0.46, 0.55, mytext)


	#plt.clim(0,6)
	#ax.set_aspect('equal')

	plt.xlabel(r'$x/L$', labelpad=1)
	plt.ylabel(r'$y/L$', labelpad=-12)
	plt.xlim(-0.5, 0.5)
	plt.ylim(-0.5, 0.5)
	plt.xticks(np.arange(-0.4, 0.4+0.1, 0.2))
	plt.yticks(np.arange(-0.4, 0.4+0.1, 0.2))

	polarization_str = '%.2f' % polarization
	mytext = '$D_r='+str(T)+', \mathcal P='+polarization_str+', t='+str(t)+'$'
#	plt.text(0.1, 0.4, mytext, color='white')
#	plt.text(0.1, 0.3, mytext, color='white')
	plt.title(mytext)

	#plt.contour(x,y, value, (0.0,0.0001), colors='k', hold='on')
	address = filename1.replace('.dat', '')
	address=address+".png"
	plt.savefig(address, transparent=True, bbox_inches='tight', pad_inches=0.01,dpi=58)
	plt.clf()
	plt.cla()
	plt.close()

arg_list = sys.argv[1:]
for i in arg_list:
	filename = i
	Plot_File(filename)


