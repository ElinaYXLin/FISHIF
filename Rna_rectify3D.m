function [Rectify]=Rna_rectify3D(Rectify_path,Rna_path,save_path,XYR,ZR)
% Rna_path='Y:\TestForEyedropper\yxt\gougouwo\hcr in chip  hb 5-7 b3 org 04_new_raw.xls';
% Rectify_path='Y:\TestForEyedropper\yxt\hcr in chip  hb 5-7 b3 org 04_new\allstack\Rectify.xml';
Rectify=readstruct(Rectify_path);
data_list = readmatrix(Rna_path);
%% xy center
data_local=data_list(:,[6 7 8]);
data_local_center=data_local;
data_local_center(:,[1 2])=data_local_center(:,[1 2])-Rectify.size/2;
%% rotate matrix
R_angle=Rectify.rotate;
Radian=(R_angle/180)*pi;
R_matrix=[...
    cos(Radian) -sin(Radian) 0;...
    sin(Radian) cos(Radian) 0;...
    0 0 1];
%% rotate
data_local_center=(R_matrix*data_local_center')';
%% 3d Rotate
RadianAxisXY=(Rectify.RotateAxisXY/180)*pi;
Zshift=Rectify.z;
if Rectify.RotateImg==1
    RadianAxisXY=-RadianAxisXY;
    Zshift=-Zshift;
end
R_matrixY=[...
    cos(RadianAxisXY(2)) 0 ZR/XYR*sin(RadianAxisXY(2));...
    0 1 0;...
    XYR/ZR*-sin(RadianAxisXY(2)) 0 cos(RadianAxisXY(2))];
R_matrixX=[...
    1 0 0;...
    0 cos(RadianAxisXY(1)) ZR/XYR*-sin(RadianAxisXY(1));...
    0 XYR/ZR*sin(RadianAxisXY(1)) cos(RadianAxisXY(1))];%Resolution fix
data_local_center=(R_matrixY*data_local_center')';
data_local_center=(R_matrixX*data_local_center')';
%% panning & Decentralisation
data_local_center(:,[1 2])=data_local_center(:,[1 2])+Rectify.xy+Rectify.size/2;
data_local_center(:,3)=data_local_center(:,3)+Zshift;
data_list(:,[6 7 8])=data_local_center;
% in figure
% data_reg=data_list(:,6)<=Rectify.size(1)&data_list(:,7)<=Rectify.size(2)...
%     &data_list(:,6)>=0&data_list(:,7)>=0;
% data_list=data_list(data_reg,:);
% out figure
data_list=data_list(:,:);
writematrix(data_list,save_path,'WriteMode','overwrite');
end

