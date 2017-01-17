from __future__ import print_function
import sys,os, subprocess
import pdb
import numpy as np
from PIL import Image

def badpixelscore(predictionfilename, truthfilename, scale, threshold=1):
	oimg = np.asarray(Image.open(predictionfilename).convert('L'))
	timg = np.asarray(Image.open(truthfilename).convert('L'))
	xmax = len(oimg[0])
	ymax = len(oimg)
	badpixels = 0
	count = ymax*xmax
	for x in range(0, xmax):
		for y in range(0, ymax):
			o = float(oimg[y][x]) # integer b/w 0 ... 255 
			t = float(timg[y][x]) # integer b/w 0 ... 255
			if (abs(o-t)>(threshold*scale)):
				badpixels += 1
	return float(badpixels)/count

def error_evaluator(hyperparameter_list):
	print('Running evaluator mode.')
	if len(sys.argv)>2:
		instancesFile = sys.argv[1]
		outputFile = sys.argv[2]
	else:
		print('USAGE: python ' + sys.argv[0] + ' <instancesFile.txt> <outputFile.csv>' )
		exit()
	baseDirectory = '../working_data/'
	with open(instancesFile, 'rb') as f:
		with open(outputFile, 'wb') as g:
			# Header for g.
			wr = 'name,orig,ae'
			for hyperparameters in hyperparameter_list:
				nranks = str(hyperparameters[0])
				ncbp = str(hyperparameters[1])
				wr += ',sym_' + nranks + '_' + ncbp
			g.write(wr+'\n')
			
			# Write g itself
			for line in f:
				s = []
				toks = line.rstrip('\n').rstrip('\r').split(',')
				name,nlabels = toks[0], toks[1]
				s.append(name)
				directory = baseDirectory + name + '/'
				print('Starting ' + name)
				truth= directory + 'disp1.png'
				orig_disp= directory + 'disp1_orig.png'
				s.append(str(badpixelscore(orig_disp, truth, 255.0/float(nlabels))))
				ae_disp= directory + 'disp1_ae.png'
				s.append(str(badpixelscore(ae_disp, truth, 255.0/float(nlabels))))
				for hyperparameters in hyperparameter_list:
					nranks = str(hyperparameters[0])
					ncbp = str(hyperparameters[1])
					sym_disp = directory + 'disp1_syms_' + nranks + '_' + ncbp + '.png'
					s.append(str(badpixelscore(sym_disp, truth, 255.0/float(nlabels))))
				g.write(','.join(s) + '\n')

# @params hyperparameter_list is a list of (nranks, ncbp) pairs
def symmetries(hyperparameter_list):
	if len(sys.argv)>1:
		instancesFile = sys.argv[1]
	else:
		print('USAGE: python ' + sys.argv[0] + ' <instancesFile.txt>' )
		exit()

	baseDirectory = '../working_data2/'
	orig_exec = './original.out'
	ae_exec = './ae.out'
	sym_exec= './sym.out'
	nranks='2'
	ncbp='0'
	with open(instancesFile, 'rb') as f:
		for line in f:
			toks = line.rstrip('\n').rstrip('\r').split(',')
			name,nlabels = toks[0], toks[1]
			directory = baseDirectory + name + '/'
			# Configuration stuff
			print('Starting ' + name)
			logfilename = directory + 'raw_symlog.txt'
			leftview = directory + 'view1.png'
			rightview= directory + 'view5.png'
			timg = directory + 'disp1.png'
			orig_disp= directory + 'disp1_orig.png'
			ae_disp= directory + 'disp1_ae.png'
			with open(logfilename,'wb') as logfile:
				# Symmetries
				for hyperparameters in hyperparameter_list:
					nranks = str(hyperparameters[0])
					ncbp = str(hyperparameters[1])
					print('\n\n\n\nSTARTING WITH nranks=' + nranks + ' and ncbp=' + ncbp + '\n\n', file=logfile)
					sym_disp = directory + 'disp1_syms_' + nranks + '_' + ncbp + '.png'
					cmd = ' '.join([sym_exec] + [leftview, rightview] + [sym_disp] + [nlabels, nranks, ncbp] + [timg])
					print('$' + cmd)
					out = subprocess.Popen( cmd , shell=True, stdout=logfile)
					out.wait()
					print('\n\n\n\nENDING WITH nranks=' + nranks + ' and ncbp=' + ncbp  + '\n\n', file=logfile)

		print('Done. exitting.')

def no_symmetries():
	if len(sys.argv)>1:
		instancesFile = sys.argv[1]
	else:
		print('USAGE: python ' + sys.argv[0] + ' <instancesFile.txt>' )
		exit()

	baseDirectory = '../working_data2/'
	orig_exec = './original.out'
	ae_exec = './ae.out'
	sym_exec= './sym.out'
	nranks='2'
	ncbp='0'
	with open(instancesFile, 'rb') as f:
		for line in f:
			toks = line.rstrip('\n').rstrip('\r').split(',')
			name,nlabels = toks[0], toks[1]
			directory = baseDirectory + name + '/'
			# Configuration stuff
			print('Starting ' + name)
			logfilename = directory + 'raw_log.txt'
			leftview = directory + 'view1.png'
			rightview= directory + 'view5.png'
			timg = directory + 'disp1.png'
			orig_disp= directory + 'disp1_orig.png'
			ae_disp= directory + 'disp1_ae.png'
			with open(logfilename,'wb') as logfile:

				# Alpha Expansion
				cmd = ' '.join([ae_exec] + [leftview, rightview] + [ae_disp] +[nlabels, nranks, ncbp] + [timg])
				print('$' + cmd)
				out = subprocess.Popen( cmd , shell=True, stdout=logfile)
				out.wait()
		print('Done. exitting.')

def hybrid(instancesFile):
	baseDirectory = '../working_data2/'
	hyb_exec = './hyb.out'
	nranks='2'
	ncbp='0'
	with open(instancesFile,'rb') as f:
		for line in f:
			toks = line.rstrip('\n').rstrip('\r').split(',')
			name,nlabels = toks[0], toks[1]
			directory = baseDirectory + name + '/'
			# Configuration stuff
			print('Starting Hybrid ' + name)
			logfilename = directory + 'raw_hybrid_log.txt'
			leftview = directory + 'view1.png'
			rightview= directory + 'view5.png'
			timg = directory + 'disp1.png'
			hyb_disp= directory + 'disp1_hyb.png'
			with open(logfilename,'wb') as logfile:
				# Alpha Expansion
				cmd = ' '.join([hyb_exec] + [leftview, rightview] + [hyb_disp] +[nlabels, nranks, ncbp] + [timg])
				print('$' + cmd)
				out = subprocess.Popen( cmd , shell=True, stdout=logfile)
				out.wait()
		print('Done. exitting.')
if __name__ == '__main__':
	print('\n\n\nSet the correct directory variables!!!!\n\n\n')
	# no_symmetries()
	# symmetries([ [1,0],[2,0],[3,0],[1,1] ])
	assert(len(sys.argv)>1),"Please provide instances list."
	hybrid(sys.argv[1])
	