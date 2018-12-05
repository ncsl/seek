import os
import sys
import numpy as np
import argparse

class ElectrodeRegistration():
    def __init__(self, subjid, brainmask_filepath, ):


if __name__ == '__main__':

    subjdir = ''

    # t1 mri filepath
    t1mri_filepath = os.path.join(subjdir, 'mri', 'orig.mgz')

    # ct file path
    ct_filepath = os.path.join(subjdir, 'ct', 'postopct.mgz')

    # brainmask data
    brainmask_filepath = os.path.join(subjdir, 'mri', 'brainmask.mgz')

    # registration matrix
    coregistration_filepath = os.path.join(subjdir, 'coregister', 'ct2t1.dat')



% datapath = '/home/pierresacre/subjects/' ;
datapath = '/media/ncsl/Pierre/subjects/';

subjID   = 1 ;
% params.letters  = {'I','T','O','J','Q','B','E','U','C','F','X','P'} ;
params.num      = [10 ,10 ,10 ,10 ,10, 10, 10, 10, 10, 10, 10, 10 ] ;
params.th = 0.75
% params.D  = 1.5+2.5+2.5 %5.5%3.5


X = orig.vol       ;
X = X(1:end,1:end-1,1:end-2) ;

B = brainmask.vol ;
B = B(1:end,1:end-1,1:end-2) ;

B_in_CT = brainmask_in_CT.vol ;
B_in_CT = B_in_CT(1:end,1:end-1,1:end-2) ;

Y = CT.vol       ;
Y = Y(1:end,1:end-1,1:end-2) ;

r_MRI = 0:size(X,1)-1; %
c_MRI = 0:size(X,2)-1; %
s_MRI = 0:size(X,3)-1; %

r_CT = 0:size(Y,1)-1; %
c_CT = 0:size(Y,2)-1; %
s_CT = 0:size(Y,3)-1; %

[R_MRI,C_MRI,S_MRI] = ndgrid(r_MRI,c_MRI,s_MRI) ;

[f,v_CRS_MRI] = isosurface(c_MRI,r_MRI,s_MRI,double(B>0)) ;

% v_CRS_CT = inv( inv(orig.tkrvox2ras) * inv(R) * CT.tkrvox2ras ) * [ v_CRS_MRI.' ; ones(1,size(v_CRS_MRI,1))] ;
v_CRS_CT = inv(CT.tkrvox2ras) * R * orig.tkrvox2ras  * [ v_CRS_MRI.' ; ones(1,size(v_CRS_MRI,1))] ;
v_CRS_CT = v_CRS_CT(1:3,:).' ;

v_RAS_CT = R * orig.tkrvox2ras  * [ v_CRS_MRI.' ; ones(1,size(v_CRS_MRI,1))] ;
v_RAS_CT = v_RAS_CT(1:3,:).' ;

