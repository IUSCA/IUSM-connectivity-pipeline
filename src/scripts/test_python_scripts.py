
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

# import os
# import numpy as np
# import nibabel as nib


# DWIpath='/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/DWI'

# fname=''.join([DWIpath,'/0_DWI.nii.gz'])
# DWI=nib.load(fname) 
# DWI_vol = DWI.get_data()

# fname=''.join([DWIpath,'/EDDY/eddy_output'])
# corrDWI=nib.load(fname)
# corrDWI_vol = corrDWI.get_data()

# corrDWI_vol = corrDWI_vol - DWI_vol


# fileOut = '/N/dc2/scratch/aiavenak/testdata/10692_1_AAK/DWI/EDDY/delta_DWI.nii.gz'
# corrDWI_new = nib.Nifti1Image(corrDWI_vol.astype(np.float32),corrDWI.affine,corrDWI.header)
# nib.save(corrDWI_new,fileOut)

# import os.path
# import numpy as np
# import nibabel as nib


# parcpath='/N/dc2/projects/connectivitypipeline/example_for_matt/matlab_outputs/T1/T1_GM_parc_yeo17_MNI152.nii'

# head_tail = os.path.split(parcpath)

# print(head_tail[0])
# print(head_tail[1])

# fileSubcort = ''.join([head_tail[0],'/T1_subcort_seg.nii'])
# print(fileSubcort)

# parc = nib.load(parcpath)
# parc_vol = parc.get_data()
# print(parc_vol.shape)
# MaxID = np.max(parc_vol)
# print(MaxID)

# subcort = nib.load(fileSubcort)
# subcort_vol = subcort.get_data()
# ind = np.argwhere(subcort_vol == 16)
# print(ind)

# subcort_vol[subcort_vol == 16] = 0
# ind = np.argwhere(subcort_vol == 16)
# print(ind)

# ids = np.unique(subcort_vol)
# print(ids)

# for s in range(0,len(ids)):
#     print(ids[s]) 
#     if ids[s] > 0:
#         subcort_vol[subcort_vol == ids[s]] = MaxID + s


# parc_vol[subcort_vol > 0] = 0
# parc_vol = np.squeeze(parc_vol)+subcort_vol

# fileOut = '/N/dc2/projects/connectivitypipeline/example_for_andrea/SUBJECTS/10692_1/T1/T1_GM_parc_yeo17_MNI152_test.nii'
# parc_vol_new = nib.Nifti1Image(parc_vol.astype(np.float32),parc.affine,parc.header)
# nib.save(parc_vol_new,fileOut)

#source activate /N/u/aiavenak/Carbonate/miniconda3/envs/scikit-env

import os.path
import numpy as np
import nibabel as nib
from skimage import measure 


EPIpath='/N/dc2/projects/connectivitypipeline/example_for_andrea/SUBJECTS/10692_1/EPI/'

fileIn="rT1_GM_parc_yeo7.nii.gz"                        
fileOut="rT1_GM_parc_yeo7_clean.nii.gz" 
thr=0.2


fileIn=''.join([EPIpath,fileIn])
v=nib.load(fileIn)  
v_vol=v.get_fdata()
print(v_vol.shape)
N = int(np.max(v_vol))
print("N = ",N)
vol_clean = np.zeros(v_vol.shape)
print("vol_clean shape ",vol_clean.shape)

for i in range(1,2):
    print(i)
    vi = v_vol == i
    vi = vi.astype(bool).astype(int)
    #print(vi)
    print("number of non-zero elements",np.count_nonzero(vi))
    clusters = measure.label(vi,neighbors=8,return_num=True)
    print("number of clusters ",clusters[1])
    for j in range(1,clusters[1]+1):
        vj = np.count_nonzero(clusters[0] == j)
        print("label ",j, "num elements ",vj)
        if vj > 8:
            print("nonzero elements in vol_clean :",np.count_nonzero(vol_clean))
            vol_clean[clusters[0] == j] = i
            print("nonzero elements in vol_clean :",np.count_nonzero(vol_clean))



fileOut=''.join([EPIpath,fileOut])

v_vol_new = nib.Nifti1Image(v_vol.astype(np.float32),v.affine,v.header)
nib.save(v_vol_new,fileOut)        



#source deactivate