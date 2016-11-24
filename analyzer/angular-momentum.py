import csv
import sys
import math


filename = sys.argv[1]

print filename
f = file(filename,'r')
data = []
n=0
for l in f:
	if (n != 0):
	    data.append(map(float,l.split()))
	n = 1
f.close()

N = math.sqrt(1+4*len(data))/2 - 1.5
print N
L = 64
M=0.0
v0=0.3

for d in data:
	if d:
		M += 2*v0*(d[0]*d[3] - d[1]*d[2]) / (N*N)

print M
