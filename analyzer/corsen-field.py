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


#import csv
#f = file("test.dat",'r')
#data = []
#for l in f:
#    data.append(map(float,l.split()))

#new_data = [[0.0 for i in xrange(5)] for j in xrange(32*32)]

#for x in xrange(32):
#	for y in xrange(32):
#		for i in xrange(5):
#			new_data[32*x + y][i] = (data[129*4*x + 4*y][i] + data[129*4*x + 4*y + 1][i] + data[129*4*x + 4*y + 2][i] + data[129*4*x + 4*y + 3][i] + data[129*(4*x + 1) + 4*y][i] + data[129*(4*x + 1) + 4*y + 1][i] + data[129*(4*x + 1) + 4*y + 2][i] + data[129*(4*x + 1) + 4*y + 3][i] + data[129*(4*x + 2) + 4*y][i] + data[129*(4*x + 2) + 4*y + 1][i] + data[129*(4*x + 2) + 4*y + 2][i] + data[129*(4*x + 2) + 4*y + 3][i] + data[129*(4*x + 3) + 4*y][i] + data[129*(4*x + 3) + 4*y + 1][i] + data[129*(4*x + 3) + 4*y + 2][i] + data[129*(4*x + 3) + 4*y + 3][i]) / 16.0

#for x in xrange(32):
#	for i in xrange(5):
#		new_data[32*x + 31][i] = (data[129*4*x + 124][i] + data[129*4*x + 125][i] + data[129*4*x + 126][i] + data[129*4*x + 127][i] + data[129*4*x + 128][i] + data[129*(4*x + 1) + 124][i] + data[129*(4*x + 1) + 125][i] + data[129*(4*x + 1) + 126][i] + data[129*(4*x + 1) + 127][i] + data[129*(4*x + 1) + 128][i] + data[129*(4*x + 2) + 124][i] + data[129*(4*x + 2) + 125][i] + data[129*(4*x + 2) + 126][i] + data[129*(4*x + 2) + 127][i] + data[129*(4*x + 2) + 128][i] + data[129*(4*x + 3) + 124][i] + data[129*(4*x + 3) + 125][i] + data[129*(4*x + 3) + 126][i] + data[129*(4*x + 3) + 127][i] + data[129*(4*x + 3) + 128][i] + data[129*(4*x + 4) + 124][i] + data[129*(4*x + 4) + 125][i] + data[129*(4*x + 4) + 126][i] + data[129*(4*x + 4) + 127][i] + data[129*(4*x + 4) + 128][i]) / 25.0

#outfile = open("output.dat",'w')
#mywriter = csv.writer(outfile, delimiter='\t')
#for x in xrange(32):
#	mywriter.writerows(new_data[x*32:(x+1)*32])
#	outfile.write('\n')

#outfile.close()

import csv
f = file("test.dat",'r')
data = []
for l in f:
    data.append(map(float,l.split()))

new_data = [[0.0 for i in xrange(5)] for j in xrange(25*25)]

for x in xrange(25):
	for y in xrange(25):
		for i in xrange(5):
			new_data[25*x + y][i] = (data[129*5*x + 5*y][i] + data[129*5*x + 5*y + 1][i] + data[129*5*x + 5*y + 2][i] + data[129*5*x + 5*y + 3][i] + data[129*5*x + 5*y + 4][i] + data[129*(5*x + 1) + 5*y][i] + data[129*(5*x + 1) + 5*y + 1][i] + data[129*(5*x + 1) + 5*y + 2][i] + data[129*(5*x + 1) + 5*y + 3][i] + data[129*(5*x + 1) + 5*y + 4][i] + data[129*(5*x + 2) + 5*y][i] + data[129*(5*x + 2) + 5*y + 1][i] + data[129*(5*x + 2) + 5*y + 2][i] + data[129*(5*x + 2) + 5*y + 3][i] + data[129*(5*x + 2) + 5*y + 4][i] + data[129*(5*x + 3) + 5*y][i] + data[129*(5*x + 3) + 5*y + 1][i] + data[129*(5*x + 3) + 5*y + 2][i] + data[129*(5*x + 3) + 5*y + 3][i] + data[129*(5*x + 3) + 5*y + 4][i) / 25.0

for x in xrange(25):
	for i in xrange(5):
		new_data[25*x + 31][i] = (data[129*5*x + 124][i] + data[129*5*x + 125][i] + data[129*5*x + 126][i] + data[129*5*x + 127][i] + data[129*5*x + 128][i] + data[129*(5*x + 1) + 124][i] + data[129*(5*x + 1) + 125][i] + data[129*(5*x + 1) + 126][i] + data[129*(5*x + 1) + 127][i] + data[129*(5*x + 1) + 128][i] + data[129*(5*x + 2) + 124][i] + data[129*(5*x + 2) + 125][i] + data[129*(5*x + 2) + 126][i] + data[129*(5*x + 2) + 127][i] + data[129*(5*x + 2) + 128][i] + data[129*(5*x + 3) + 124][i] + data[129*(5*x + 3) + 125][i] + data[129*(5*x + 3) + 126][i] + data[129*(5*x + 3) + 127][i] + data[129*(5*x + 3) + 128][i] + data[129*(5*x + 4) + 124][i] + data[129*(5*x + 4) + 125][i] + data[129*(5*x + 4) + 126][i] + data[129*(5*x + 4) + 127][i] + data[129*(5*x + 4) + 128][i]) / 25.0

outfile = open("output.dat",'w')
mywriter = csv.writer(outfile, delimiter='\t')
for x in xrange(25):
	mywriter.writerows(new_data[x*25:(x+1)*25])
	outfile.write('\n')

outfile.close()
