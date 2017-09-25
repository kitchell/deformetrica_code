import os
import numpy as np
import nibabel as nib
from PIL import Image

#Sizes used in the functions.
xlen = 128
ylen = 128
zlen = 128

#Generates a .nii file with a white ellipsoid in a black background.
def generateEllipsoidImage(a,b,c):
  data = np.zeros((xlen, ylen, zlen))
  for x in range(xlen):
    for y in range(ylen):
      for z in range(zlen):
        if ((x-xlen*1./2.)/a)**2 + ((y-ylen*1./2.)/b)**2 + ((z-zlen*1./2.)/c)**2<=1:
          data[x,y,z] = 1
  img = nib.Nifti1Image(data, np.eye(4))
  nib.save(img, "ellipsoidImage"+str(a)+"_"+str(b)+"_"+str(c)+".nii")

#Generates a .nii file with a white square in a black background.
def generateSquareImage(a):
     data = np.zeros((xlen, ylen, zlen))
     for x in range(xlen):
       for y in range(ylen):
         for z in range(zlen):
           if (abs(x-xlen*1./2)<=a*1./2 and abs(y-ylen*1./2)<=a*1./2 and abs(z-zlen*1./2)<=a*1./2 ):
             data[x,y,z] = 1;
     img = nib.Nifti1Image(data, np.eye(4))
     nib.save(img, "ellipsoidSquare"+str(a)+str(xlen)+".nii")

#Generates a .png file with a rectangle square in a black background.
def generateRectangleImage(a, b):
     data = np.zeros((xlen, ylen))
     for x in range(xlen):
       for y in range(ylen):
           if (abs(x-xlen*1./2)<=a and abs(y-ylen*1./2)<=b):
             data[x,y] = 256
     img = Image.fromarray(data)
     img.convert("L").save("rectangleImage"+str(a)+"_"+str(b)+".png")

#Generates a .png file with an ellipse in a black background.
def generateEllipseImage(a,b):
    data = np.zeros((xlen, ylen))
    for x in xrange(xlen):
     for y in xrange(ylen):
         if ((x-xlen*1./2.)/a)**2 + ((y-ylen*1./2.)/b)**2 <= 1:
             data[x,y] = 256
    img = Image.fromarray(data)
    img.convert("L").save("ellipseImage"+str(a)+"_"+str(b)+".png")

