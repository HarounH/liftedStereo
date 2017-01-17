from __future__ import print_function

'''
g++ hoStereo4.cpp -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -lexternal-library-qpbo-shared -o a.out
'''

import sys, os, subprocess
executable = './a.out'

# Step 1 - setup parameters
# Params - image
imagenumber = '1'
lfile = '../data/'+imagenumber+'l.png'
rfile = '../data/'+imagenumber+'r.png'
tfile = '../data/'+imagenumber+'t.png'
# ASSERT  output files aren't ready yet.

# Params - problem parameters
nLabels = '32'
goldScale= '8'
badPixelThreshold='2.0'

# Energy parameters
unaryPotentialRadius='3'
truncationThreshold = '2'
truncationPenalty='20'

# Inference parameters
nIterations='50'

# Symmetry parameters
nRanksUsed = '4' # TODO alter this and figure out performance.
nIterationsOfCBP='1'

# Output files
nosyms_ofile = '../data/reductions/' + imagenumber + 'o_ns'+ nIterations +'_' + nLabels + '_' + truncationThreshold + '_' + truncationPenalty + '.png'
syms_ofile = '../data/reductions/' + imagenumber + 'o_s'+ nIterations +'_' + nLabels + '_' + truncationThreshold + '_' + truncationPenalty + '_' + nRanksUsed + '_' + nIterationsOfCBP + '.png'
# Step 2 - Setup command and streams.
cmd = ' '.join([executable, lfile, rfile, nosyms_ofile, syms_ofile, tfile, nLabels, goldScale, badPixelThreshold, unaryPotentialRadius, truncationThreshold, truncationPenalty, nIterations, nRanksUsed, nIterationsOfCBP])
# Step 3 - call command
with open('../logs/' + imagenumber + 'error_plots' + '0.txt','wb') as f:
	out = subprocess.Popen( cmd , shell=True, stdout=f)
	out.wait()
