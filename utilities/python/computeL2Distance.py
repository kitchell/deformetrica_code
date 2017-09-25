import os
import numpy as np
import nibabel as nib
import sys
from PIL import Image


assert len(sys.argv)>=3, "Usage: computeL2Distance.py image1.nii image.nii"

imagePath1 = sys.argv[1]
imagePath2 = sys.argv[2]


#print("Reading image from files : ", imagePath1, imagePath2)
img1 = nib.load(imagePath1)
data1 = img1.get_data()
x1,y1,z1 = data1.shape
del img1
min1 = np.min(data1)
max1 = np.max(data1)
data1 = (data1-min1)/max1
#print("Image shape :", x1, y1, z1)
#print("min intensity image1 : ", min1, "max intensity image 1 : ", max1)
img2 = nib.load(imagePath2)
data2 = img2.get_data()
x2,y2,z2 = data2.shape
del img2
min2 = np.min(data2)
max2 = np.max(data2)
data2 = (data2-min2)/max2 #this is LPS
#print("min intensity image2 : ", min2, "max intensity image 2 : ", max2)

assert ((x1==x2)and(y1==y2)and(z1==z2)), "Images dimensions mismatch (x1,y1,z1,x2,y2,z2) %i %i %i vs %i %i %i" % (x1,y1,z1,x2,y2,z2)

L2 = 0
L2 = np.linalg.norm(data1-data2)
print("L2 Distance : ", L2, "(per pixel :)",L2/(x1*y1*z1))
