import sys,os
import numpy as np
from PIL import Image

def normalize(infile, outfile):
	img = np.asarray(Image.open(infile).convert('L'))
	mi = img.min()
	ma = img.max()
	mul = 254.0/(ma-mi)
	Image.fromarray(((img-mi + 1)*mul).astype(np.uint8)).save(outfile)

normalize(sys.argv[1], sys.argv[2])