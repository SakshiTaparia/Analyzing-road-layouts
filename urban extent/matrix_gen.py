import numpy as np
import scipy
from scipy.misc import *
from scipy import misc
from scipy import ndimage
import matplotlib
from matplotlib.pyplot import *
import PIL
from PIL import *
#district = 'Gurgaon_summer'
#dataset = matplotlib.pyplot.imread('delhi_bu_nbu.tif')
#dataset *= 120
#imsave('output.png', dataset)
img = np.array(Image.open('bangalore_new.png'))
f = open('new_input.txt','w')
for i in range(0,img.shape[0]):
	for j in range(0,img.shape[1]):
		f.write(str(img[i][j]))
		f.write(' ')
	f.write('\n')
img*=120
img2 = Image.fromarray(img)
img2.show()
#1858 * 2228