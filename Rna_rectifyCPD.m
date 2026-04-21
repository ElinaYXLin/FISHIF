function [Rectify]=Rna_rectifyCPD(stack_path,Rna_path,save_path)
% Rna_path='Y:\TestForEyedropper\yxt\gougouwo\hcr in chip  hb 5-7 b3 org 04_new_raw.xls';
% Rectify_path='Y:\TestForEyedropper\yxt\hcr in chip  hb 5-7 b3 org 04_new\allstack\Rectify.xml';
load([stack_path,'Transform_m.mat']);
data_list = readmatrix(Rna_path);
%% get size
stack_name=dir([stack_path,'stack*.tif']);
if size(stack_name,1)==0
    stack_name=dir([stackpath,'time0001stack*.tif']);
end
z_dim=size(stack_name,1);
stack_name=stack_name.name;
img=imread([stack_path,stack_name]);
Rectify.size=size(img(:,:,1));
%% local
data_local=data_list(:,[6 7 8]);
%% rotate matrix
data_local_rectify=(R_m*data_local'+T_m)';
data_list(:,[6 7 8])=data_local_rectify;
data_reg=data_list(:,6)<=Rectify.size(1)&data_list(:,7)<=Rectify.size(2)&data_list(:,8)<=z_dim...
    &data_list(:,6)>=0&data_list(:,7)>=0&data_list(:,8)>=0;
data_list=data_list(data_reg,:);
writematrix(data_list,save_path,'WriteMode','overwrite');
end