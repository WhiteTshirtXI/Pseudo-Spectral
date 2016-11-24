#f = file("test.dat",'r')
#data = []
#for l in f:
#    data.append(map(float,l.split()))

#new_data = [[0.0 for i in xrange(5)] for j in xrange(64*64)]

#for x in xrange(64):
#	for y in xrange(64):
#		for i in xrange(5):
#			new_data[64*x + y][i] = (data[129*2*x + 2*y][i] + data[129*2*x + 2*y + 1][i] + data[129*(2*x + 1) + 2*y][i] + data[129*(2*x + 1) + 2*y + 1][i]) / 4.0

#for x in xrange(64):
#	for i in xrange(5):
#		new_data[64*x + 63][i] = (data[129*2*x + 126][i] + data[129*2*x + 127][i] + data[129*2*x + 128][i] + data[129*(2*x + 1) + 126][i] + data[129*(2*x + 1) + 127][i] + data[129*(2*x + 1) + 128][i]) / 6.0

#outfile = open("output.dat",'w')
#mywriter = csv.writer(outfile, delimiter='\t')
#for x in xrange(64):
#	mywriter.writerows(new_data[x*64:(x+1)*64])
#	outfile.write('\n')

#outfile.close()


import csv
import math

f = file("test.dat",'r')
data = []
for l in f:
    data.append(map(float,l.split()))

f.close()

relation = [0 for i in xrange(129*128)]

for x in xrange(129):
	for y in xrange(129):
		if (data[x*128+y][4] <= 0):
			data[x*128+y][4] = 0.00001

for x in xrange(129):
	for y in xrange(128):
		relation[x*128+y] = math.log(data[(128-x)*129+y+1][4] / data[(128-x)*129+y][4])

xc = []
yc = []

for x in xrange(1,129):
	for y in xrange(1,128):
		if (math.fabs(relation[x*128+y]) > 0.1):
			xc.append(float(x)/129.0)
			yc.append(float(y)/129.0)
			break

for x in xrange(len(xc)):
		print xc[x], yc[x]

