import sys,os
import matplotlib.pyplot as plt
import numpy as np

bases = ['art','baby1', 'baby2','baby3']
aegm_files = ['aegm-' +base + '.csv' for base in bases]
sym10_files= ['sym-1-0-' +base + '.csv' for base in bases]
sym11_files= ['sym-1-1-' +base + '.csv' for base in bases]
sym20_files= ['sym-2-0-' +base + '.csv' for base in bases]
sym21_files= ['sym-2-1-' +base + '.csv' for base in bases]
sym30_files= ['sym-3-0-' +base + '.csv' for base in bases]

aegm_iter=[]
aegm_energy=[]
aegm_time=[]

sym10_iter=[]
sym10_energy=[]
sym10_time=[]

sym11_iter=[]
sym11_energy=[]
sym11_time=[]

sym20_iter=[]
sym20_energy=[]
sym20_time=[]

sym21_iter=[]
sym21_energy=[]
sym21_time=[]

sym30_iter=[]
sym30_energy=[]
sym30_time=[]

mode = sys.argv[1] # Either 
if mode=='iter':
	xIdx = 0
elif mode=='time':
	xIdx = 2
else:
	print('USAGE: python ' + sys.argv[0] + '<iter or time>')
	exit()

# Update syms stuff.
def fetchData(filenames, iterData, energyData, timeData):
	iters = []
	energies=[]
	times=[]
	for file in filenames:
		itr = []
		energy=[]
		time=[]
		with open(file) as f:
			mi = float(f.readline().split(',')[1])
			for line in f:
				toks = line.split(',')
				itr.append(int(toks[0]))
				energy.append(float(toks[1])/mi)
				time.append(float(toks[2]))
		iters.append(itr)
		energies.append(energy)
		times.append(time)
	# Time to accumulate the data into iterData, energyData, timeData
	nSamples = len(filenames)
	for i in range(0, max([len(t) for t in iters])):
		iterData.append(0.0)
		energyData.append(0.0)
		timeData.append(0.0)
		cnt = 0
		for j in range(0, nSamples):
			if len(iters[j]) > i:
				cnt += 1
				iterData[-1] += iters[j][i]
				energyData[-1] += energies[j][i]
				timeData[-1] += times[j][i]
		if cnt==0:
			continue
		iterData[-1] = iterData[-1]/float(cnt)
		energyData[-1] = energyData[-1]/float(cnt)
		timeData[-1] = timeData[-1]/float(cnt)

def addPlot(mode, itr, energy, time, label):
	if mode=='iter':
		plt.plot(itr, energy, label=label)
	elif mode == 'time':
		plt.plot(time, energy, label=label)
	else:
		print "unknown mode for addPlot. exitting"
		exit()

fetchData(aegm_files, aegm_iter, aegm_energy, aegm_time)
fetchData(sym10_files, sym10_iter, sym10_energy, sym10_time)
fetchData(sym11_files, sym11_iter, sym11_energy, sym11_time)
fetchData(sym20_files, sym20_iter, sym20_energy, sym20_time)
fetchData(sym21_files, sym21_iter, sym21_energy, sym21_time)
fetchData(sym30_files, sym30_iter, sym30_energy, sym30_time)

addPlot(mode, aegm_iter, aegm_energy, aegm_time, 'aegm')
addPlot(mode, sym10_iter, sym10_energy, sym10_time, 'sym10')
addPlot(mode, sym20_iter, sym20_energy, sym20_time, 'sym20')
addPlot(mode, sym11_iter, sym11_energy, sym11_time, 'sym11')
addPlot(mode, sym21_iter, sym21_energy, sym21_time, 'sym21')
addPlot(mode, sym30_iter, sym30_energy, sym30_time, 'sym30')

plt.title('energy/min vs ' + mode)
plt.xlabel(mode)
plt.ylabel('energy/min')
plt.legend()
plt.show()