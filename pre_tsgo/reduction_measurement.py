from __future__ import print_function

'''
g++ hoStereo1.cpp -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -lexternal-library-qpbo-shared -o a.out
g++ hoStereo2.cpp -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -lexternal-library-qpbo-shared -o b.out
g++ hoStereo3.cpp -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -lexternal-library-qpbo-shared -o c.out
'''

import sys, os, subprocess

# Trying to measure the reduction of model using different nRanksToLook at and stuff.


map_executable = './c.out'
nosyms_executable='./a.out'
# symsExecutable = 
mode='single'


nFineness='2'
nMaxIterations='1'

imagenumber = '1'
nLabels='32'
lfile = '../data/'+imagenumber+'l.png'
rfile = '../data/'+imagenumber+'r.png'
tfile = '../data/'+imagenumber+'t.png'
sfile = '../data/reductions/ignore.png'
truncationThreshold='1'
# lmbda= '20'
# pnLambda='20'
lmbda = str(20)
pnLambda=str(10)
upr='1'
nSteps='60'

print('\n\n\n\nDONT FORGET TO RECOMPILE!!!!\n\n\n\n')

# Symmetries
ofile = '../data/reductions/'+imagenumber+'o_qpbo_'+nLabels+'_'+lmbda+'_map_' + nFineness + '_' + nMaxIterations + '.png'
sfile = '../data/reductions/syms/'+imagenumber+'s_qpbo_'+nLabels+'_'+lmbda+'_map_' + nFineness + '_' + nMaxIterations + '.png'
cmd = ' '.join([map_executable, mode, lfile, rfile, ofile, nLabels, truncationThreshold, lmbda, pnLambda, upr, nSteps, tfile, nFineness, nMaxIterations, sfile])
print(cmd)
with open('../logs/get_groups.txt','wb') as f:
	out = subprocess.Popen( cmd , shell=True, stdout=f)

# No symmetries
# ofile = '../data/reductions/'+imagenumber+'o_qpbo_'+nLabels+'_'+lmbda+'_syms.png'
# sfile = '../data/reductions/syms/'+imagenumber+'s_qpbo_'+nLabels+'_'+lmbda+'_syms.png'
# cmd = ' '.join([nosyms_executable, mode, lfile, rfile, ofile, nLabels, truncationThreshold, lmbda, pnLambda, upr, nSteps, tfile])
# out = subprocess.Popen(cmd, shell=True)
# print(cmd)

out.wait()
exit()

with open('../logs/map_model_reduction_on'+imagenumber+'_'+nLabels+'_'+lmbda+'_qpbo.txt', 'wb') as f:
	with open('../logs/nosyms_accuracies_on'+imagenumber+'_'+nLabels+'_'+lmbda+'_qpbo.txt', 'wb') as g:	
		for fineness in range(1, 3):
			for nCBPIterations in range(0,2):
				nFineness = str(fineness)
				nMaxIterations = str(nCBPIterations)
				ofile = '../data/reductions/'+imagenumber+'o_qpbo_'+nLabels+'_'+lmbda+'_map_' + nFineness + '_' + nMaxIterations + '.png'
				sfile = '../data/reductions/syms/'+imagenumber+'s_qpbo_'+nLabels+'_'+lmbda+'_map_' + nFineness + '_' + nMaxIterations + '.png'
				cmd = ' '.join([map_executable, mode, lfile, rfile, ofile, nLabels, truncationThreshold, lmbda, pnLambda, upr, nSteps, tfile, nFineness, nMaxIterations, sfile])
				print('Now:' + nFineness + ' ' + nMaxIterations + ' starting:')
				print(cmd)
				out = subprocess.Popen( cmd , shell=True, stdout=f)
				out.wait()
		ofile = '../data/reductions/'+imagenumber+'o_qpbo_'+nLabels+'_'+lmbda+'_nosyms.png'
		sfile = '../data/reductions/syms/'+imagenumber+'s_qpbo_'+nLabels+'_'+lmbda+'_nosyms.png'
		cmd = ' '.join([nosyms_executable, mode, lfile, rfile, ofile, nLabels, truncationThreshold, lmbda, pnLambda, upr, nSteps, tfile])
		print('Now: nosyms starting:')
		print(cmd)
		out = subprocess.Popen( cmd , shell=True, stdout=g)
		out.wait()