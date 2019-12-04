
#!/bin/python

# ------------------------------------------------------------
##### get_private_tags

# import pydicom
# import sys

# path2dicom=str(sys.argv[1])
# print(path2dicom)

# dicomHeader=pydicom.read_file(path2dicom)
# matrix=dicomHeader[0x0051100b].value
# dim1=matrix.split('*')
# dim1=int(dim1[1])

# # accelerator factor
# try: 
#     AcqStr=dicomHeader[0x0051100b].value
#     Pstring=AcqStr.split('p')
#     AccF=int(Pstring[0])
# except:
#     AccF=1
    
# # bandwidth per pixel phase encoding (hz)
# bppe=int(dicomHeader[0x00191028].value)

# # Effective Echo Spacing (s)
# ees=1/(bppe*dim1)

# #actual number of phase encoding lines
# anofel=dim1/AccF

# # Total Readout Time (s)
# RT=(anofel-1)*ees
# #print('RT ',RT)
# sys.stdout.write(str(RT))

# # ------------------------------------------------------------
# ##### get_largest_clusters_only

# import nibabel as nib
# import numpy as np

# vol=nib.load('/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/rT1_GM_parc_shen_278.nii.gz')

# print(vol.shape)
# # print(vol)

# vol_data=vol.get_fdata()

# print(np.unique(vol_data))

# # ------------------------------------------------------------
# ##### read a mask

# import nibabel as nib
# import numpy as np

# resting=nib.load('/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz')
# resting_vol=resting.get_data()
# print(resting.shape)
# [sizeX,sizeY,sizeZ,numTimePoints] = resting_vol.shape

# # read brain mask
# volBrain=nib.load('/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/rT1_brain_mask.nii.gz')
# print(volBrain.shape)
# volBrain_vol=volBrain.get_data()

# volRef=nib.load('/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/2_epi_meanvol_mask.nii.gz')
# print(volRef.shape)
# volRef_vol=volRef.get_data()

# volBrain_vol = (volBrain_vol > 0) & (volRef_vol != 0)
# print(volBrain_vol.shape)
# ones_vol=np.sum(volBrain_vol)
# print(ones_vol)

# fileOut = '/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/rT1_brain_mask_FC.nii.gz'
# volBrain = nib.Nifti1Image(volBrain_vol.astype(np.float32),volBrain.affine())
# nib.save(volBrain,fileOut)

# # ------------------------------------------------------------
# ##### demean and detrend

# import nibabel as nib
# import numpy as np
# from scipy import signal

# resting=nib.load('/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz')
# resting_vol=resting.get_data()
# print(resting.shape)
# [sizeX,sizeY,sizeZ,numTimePoints] = resting_vol.shape

# # read brain mask
# volBrain=nib.load('/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/rT1_brain_mask.nii.gz')
# print(volBrain.shape)
# volBrain_vol=volBrain.get_data()

# for i in range(0,sizeX-1):
#     for j in range(0,sizeY-1):
#         for k in range(0,sizeZ-1):
#             if volBrain_vol[i,j,k] > 0:
#                 TSvoxel = resting_vol[i,j,k,:].reshape(numTimePoints,1)
#                 mean = np.mean(TSvoxel)
#                 TSvoxel_detrended = signal.detrend(TSvoxel,type='linear')
#                 resting_vol[i,j,k,:] = TSvoxel_detrended.reshape(1,1,1,numTimePoints)
#     if i % 25 == 0:
#         print(i/sizeX)  ## change this to percentage progress 

# print(TSvoxel_detrended)
# fileOut2 = '/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/6_epi.nii.gz'
# resting_detrended = nib.Nifti1Image(resting_vol.astype(np.float32),resting.affine,resting.header)
# nib.save(resting_detrended,fileOut2)


# fileOut = '/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/EPI/rT1_brain_mask_FC.nii.gz'
# volBrain = nib.Nifti1Image(volBrain_vol.astype(np.float32),volBrain.affine)
# nib.save(volBrain,fileOut)

# # ------------------------------------------------------------
# ##### DWI_A

# import os
# from dipy.io import read_bvals_bvecs
# import nibabel as nib
# import numpy as np

# # p=os.environ['path']
# p='/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/DWI'

# pbval=''.join([p,'/0_DWI.bval'])
# pbvec=''.join([p,'/0_DWI.bvec'])

# bvals, bvecs = read_bvals_bvecs(pbval,pbvec)
# # print("bvals size", bvals.shape)
# # print("bvecs size", bvecs.shape)

# if bvals.shape[0] > 1:
#     # vector is vertical, needs to be transposed
#     bvals = bvals.reshape((1,bvals.size)) 
#     # print("bvals size", bvals.shape)

# if bvecs.shape[0] > 3:
#      # vector is vertical, needs to be transposed
#     bvecs = bvecs.T 
#     # print("bvecs size", bvecs.shape)

# DWIp=''.join([p,'/0_DWI.nii.gz'])
# DWI=nib.load(DWIp)  
# # DWI_vol=DWI.get_data()
# # print(DWI.shape)

# pbvalt=''.join([p,'/0_DWIt.bval'])
# pbvect=''.join([p,'/0_DWIt.bvec'])

# # print('bvals.shape[1] ',bvals.shape[1])
# # print('bvecs.shape[1] ',bvecs.shape[1])
# # print('DWI.shape[3] ',DWI.shape[3])

# if bvals.shape[1] == DWI.shape[3] and bvecs.shape[1] == DWI.shape[3]:
#     np.savetxt(pbvalt,bvals,delimiter='\t',fmt='%u')
#     np.savetxt(pbvect,bvecs,delimiter='\t',fmt='%f')
#     print('1')
# else:
#     print('0')

# # ------------------------------------------------------------
# ##### DWI_A - topoup

# import os
# import numpy as np


# def is_empty(any_struct):
#     if any_struct:
#         return False
#     else:
#         return True 

# DWIpath='/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/DWI'

# pbval=''.join([DWIpath,'/0_DWI.bval'])
# bval = np.loadtxt(pbval)
# print(bval)

# B0_index = np.where(bval<=1)
# print(B0_index)

# if is_empty(B0_index):    
#     #print("No B0 volumes identified. Check quality of 0_DWI.bval") 
#     print(0)
# else:     
#     ff = open("test.txt","w+")
#     for i in np.nditer(B0_index):
#         fn = "/AP_b0_%d.nii.gz" % i
#         fileOut = ''.join([DWIpath,fn])
#         ff.write("%s \n" % fileOut)
#         print(fileOut)
#     ff.close()
#     print(1)

#     # it = np.nditer(B0_index,flags=['f_index'])
#     # while not it.finished:
#     #     fn = "/AP_b0_%d.nii.gz" % it[0]
#     #     fileOut = ''.join([DWIpath,fn])
#     #     print(fileOut)
#     #     #print(it[0])
#     #     #print(it.index)
#     #     it.iternext()      

# # # ------------------------------------------------------------
# # ##### DWI_A - EDDY

# import os
# import numpy as np
# import nibabel as nib

# def read_integers(filename):
#     with open(filename) as f:
#         return map(int, f)    

# def is_empty(any_struct):
#     if any_struct:
#         return False
#     else:
#         return True 

# DWIpath='/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/DWI'

# # read in DWI data and find number of volumes
# fname=''.join([DWIpath,'/0_DWI.nii.gz'])
# DWI=nib.load(fname)  
# ss=DWI.shape
# numVols=ss[3]
# #print("numvols = ",numVols)

# b0file = ''.join([DWIpath,"/b0file.txt"])

# ff = open(b0file,"r")
# ffl = ff.readlines()

# Index=np.ones((numVols,1),dtype=np.int64)

# for i in range(0,len(ffl)):
#     ii = int(ffl[i]) 
#     if ii != 1:  
#         #  for every subsequent B0 the volume index increases. This provides temporal information about location of B0 volumes
#         Index[ii:]=i+1


# # save to file
# fname=''.join([DWIpath,'/EDDY/index.txt'])
# np.savetxt(fname,Index, fmt='%s')

# ff.close()
# print(1)
# # ------------------------------------------------------------
# ##### DWI_A - EDDY- last section

import os
import numpy as np
import nibabel as nib


DWIpath='/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/DWI'

fname=''.join([DWIpath,'/0_DWI.nii.gz'])
DWI=nib.load(fname) 
DWI_vol = DWI.get_data()

fname=''.join([DWIpath,'/EDDY/eddy_output'])
corrDWI=nib.load(fname)
corrDWI_vol = corrDWI.get_data()

corrDWI_vol = corrDWI_vol - DWI_vol


fileOut = '/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/DWI/EDDY/delta_DWI.nii.gz'
corrDWI_new = nib.Nifti1Image(corrDWI_vol.astype(np.float32),corrDWI.affine,corrDWI.header)
nib.save(corrDWI_new,fileOut)