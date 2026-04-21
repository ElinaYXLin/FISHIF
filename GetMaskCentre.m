function [Centroid_nucleus]=GetMaskCentre(DAPI_channel,varargin)
%% Program Intro
% a program to get nucleus centres by mask.mat
% GetMaskCentre(DAPI_channel,{maskpath,stackpath})
%% Get the path of mask file 
if isempty(varargin)
    [~,maskpath,~]=uigetfile({'*.mat'},'Select A Matlab Data');
    [~,stackpath,~]=uigetfile('.tif','Select A Stack');
else
    maskpath=varargin{1}{1};
    stackpath=varargin{1}{2};
end
%% 
mask_nmae  = dir([maskpath,'*mask.mat']);
time_dim = size(dir([stackpath,'*stack01.tif']),1);
for time_index=1:time_dim
    try
        load([maskpath,'time',num2str(time_index,'%04u'),'mask.mat']);
    end
    stack_name=dir([stackpath,'time',num2str(time_index,'%04u'),'stack*.tif']);
    if size(stack_name,1)==0 && time_dim==1
        stack_name=dir([stackpath,'stack*.tif']);
        load([maskpath,'mask.mat']);
    end
    stack_name={stack_name.name};
    nucleus_num=double(max(mask_stack(:)));
    stack_num=size(mask_stack,3);
    MeanIntensity_layers=nan(nucleus_num,stack_num);
    Centroid_layers=cell(nucleus_num,stack_num);
    for I_layer = 1:size(mask_stack,3)
        img=imread([stackpath,stack_name{I_layer}]);
        sg_prop = regionprops(double(mask_stack(:,:,I_layer)),double(img(:,:,DAPI_channel)),'MeanIntensity','Centroid');
        MeanIntensity_layers(1:size(sg_prop,1),I_layer)=[sg_prop.MeanIntensity]';
        Centroid_layers(1:size(sg_prop,1),I_layer)={sg_prop.Centroid}'; %[x y]
    end
    
    [~,max_index]=max(MeanIntensity_layers,[],2);
    Centroid_layers_focus=Centroid_layers(:,max_index');
    Centroid_layers_focus=Centroid_layers_focus(1:(nucleus_num+1):(nucleus_num*nucleus_num))';
    Centroid_layers_focus=cell2mat(Centroid_layers_focus);
%     Centroid_nucleus_index=[Centroid_layers_focus,max_index];
%     Centroid_nucleus=[size(mask_stack,1)-Centroid_layers_focus(:,2),Centroid_layers_focus(:,1),max_index];%[row col z]
    Centroid_nucleus=[Centroid_layers_focus(:,2),Centroid_layers_focus(:,1),max_index];%[row col z]
end
end