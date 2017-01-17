from __future__ import print_function

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
output= sys.argv[1]
truth = sys.argv[2]
nlabels= int(sys.argv[3])
scale = int(sys.argv[4])
threshold = 1
# Iterate over all pixels.

oimg = np.asarray(Image.open(output).convert('L'))
timg = np.asarray(Image.open(truth).convert('L'))

xmax = len(oimg[0])
ymax = len(oimg)

badpixels = 0
count = ymax*xmax
for x in range(0, xmax):
	for y in range(0, ymax):
		o = oimg[y][x] # Its a triplet....but what's the disparity?
		od = nlabels*float(o)/255 # Thats the disparity.
		t = timg[y][x]
		td = float(t)/scale

		if abs(od-td)>threshold:
			badpixels += 1
print('error=' + str(float(badpixels*100)/count) + '%')