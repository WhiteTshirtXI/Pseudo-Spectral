import csv
import sys
import math


def Find_Polarity(filename):
	f = file(filename,'r')
	data = []
	n=0
	for l in f:
		if (n != 0):
		    data.append(map(float,l.split()))
		n = 1
	f.close()

	N = round(math.sqrt(len(data) - 3.0/4.0) - 0.5)
#	print N
	Px=0.0
	Py=0.0

	for d in data:
		if d:
			Px += d[2] / (N*N)
			Py += d[3] / (N*N)

	P2 = Px*Px + Py*Py
	P = math.sqrt(P2)

	start = filename.find('T=') + 2
	end = filename.find('.dat', start)
	print filename[start:end],P

arg_list = sys.argv[1:]
for i in arg_list:
	filename = i
	Find_Polarity(filename)
