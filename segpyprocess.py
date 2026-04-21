# Jianing Zhu
# exploit time: 2021/9/13 09:21
import numpy as np
import scipy.io as scio
import os
import pandas as pd
def mask_stack(filepath,savepath):
    '''

    :param filepath: str,input _seg.npy filefolder
    :return:
    '''
    filename_npy = []
    file_npy = []
    os.chdir(filepath)
    file_chdir = os.getcwd()


    for root,dirs,files in os.walk(file_chdir):
        for file in files:
            if os.path.splitext(file)[-1] == '.npy':

                 dat = np.load(file, allow_pickle=True).item()
                 scio.savemat(str(savepath+'\mask_'+os.path.splitext(file)[0]+'.mat'), {'masks': dat['masks']})
                #scio.savemat(str(savepath + '\est_diam_' + os.path.splitext(file)[0] + '.mat'), {'est_diam': dat['est_diam']})
                 filename_npy.append(file)
                 file_npy.append(np.load(file, allow_pickle=True))
    maskdata = file_npy
    return maskdata
#
a = mask_stack('E:\\data\\wjy_data\\10222021\\10212021_PE_hb1_tmr_hb7_647_yell_HCR488_60X_2021_10_21__21_22_43_005_new\\raw2input','E:\\data\\wjy_data\\masks\\10222021\\10212021_PE_hb1_tmr_hb7_647_yell_HCR488_60X_2021_10_21__21_22_43_005_new')







    # for i in range():
    #     data = np.load('filepath',allow_pickle=True).item()



# dat = np.load('time0003stack02_seg.npy', allow_pickle=True).item()


# data_address = 'E:\data\mask1.mat'
# scio.savemat(data_address,{'masks':dat['masks']})