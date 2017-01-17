import sys,os
import matplotlib.pyplot as plt
import numpy as np

csvnames = ['aegm-art.csv', 'sym-1-0-art.csv','sym-1-1-art.csv','sym-2-0-art.csv'] # Fill this in

mode = sys.argv[1] # Either 
if mode=='iter':
	xIdx = 0
elif mode=='time':
	xIdx = 2
else:
	print('USAGE: python ' + sys.argv[0] + '<iter or time>')
	exit()

for file in csvnames:
	with open(file, 'rb') as f:
		topline = f.readline().split(',')
		mi = float(topline[1])
		x = []
		energy = []
		for line in f:
			toks = line.split(',')
			x.append(float(toks[ xIdx ]))
			energy.append(float(toks[1])/mi)
			# time.append(float(toks[2]))
		plt.plot(x, energy, label=(file.split('.')[0]))

plt.title('energy/min vs ' + mode)
plt.xlabel(mode)
plt.ylabel('energy/min')
plt.legend()		
plt.show()