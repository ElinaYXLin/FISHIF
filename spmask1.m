function [mask_out,varargout] = spmask1(imstack,Nr,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% A function for spot recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mask_out (Output): spot candidates;
%% varargout (Output): mask_out2D: 2D maxima;
%% imstack (input): input image stack;
%% Nr (input): spot spread half width (on z direction);
%% varargin (input): (1) im_max: peak spot intensity threshold (default = 0);
%%                   (2) Nr3: 3D maxima extension range for 2D maxima recognition (default = 2);
%%                   (3) z00: boundary layers that do not need spot detection (default = 1);
%%                   (4) Nr_or: whether to use "OR" for 3D spot spread half width (Nr) check (default = false);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) && ~isempty(varargin{1})
    im_max = varargin{1};
else
    im_max = 0;
end

if length(varargin) > 1 && ~isempty(varargin{2})   %%% 3D maxima extension range
    Nr3 = ceil(varargin{2});
else
    Nr3 = 2;
end

if length(varargin) > 2 && ~isempty(varargin{3})   %%% boundary layers that do not need spot detection
    z00 = varargin{3};
else
    z00 = 1;
end

if length(varargin) > 3 && ~isempty(varargin{4})   %%% boundary layers that do not need spot detection
    Nr_or = varargin{4};
else
    Nr_or = false;
end

Ixy = 3;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


max3 = imregionalmax(imstack);   %%% 3D local maxima
max3(:,:,[1:z00,end+1-z00:end]) = false;

max30 = false(size(max3));   %%% 2D local maxima of 3D local maxima residues
max2 = false(size(max3));   %%% 2D local maxima
max20 = max2;
N_layer = size(imstack,3);
for I_layer = 1:N_layer
    max20(:,:,I_layer) = imregionalmax(imstack(:,:,I_layer));
    max2(:,:,I_layer) = imdilate(max20(:,:,I_layer),strel('square',Ixy));
end

Nr = ceil(Nr);
Ar = -Nr:Nr;
zmin = max(1,(1+Ar));
zmax = min(N_layer,(N_layer+Ar));
if Nr_or
    for Ir = (length(Ar)+1)/2:length(Ar)
        Ir2 = length(Ar)-Ir+1;
        temp = max3;
        temp(:,:,zmin(Ir2):zmax(Ir2)) = temp(:,:,zmin(Ir2):zmax(Ir2)) & max2(:,:,zmin(end-Ir2+1):zmax(end-Ir2+1));
        max3(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
        max3 = max3 | temp;
    end
else
    for Ir = 1:length(Ar)
        max3(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
    end
end

cum_max2 = (cumsum(~max2,3)+1).*max2;
num_max2 = zeros([size(max2(:,:,1)),max(cum_max2(:))]);
for ii = 1:max(cum_max2(:))
    num_max2(:,:,ii) = sum(cum_max2 == ii,3);
end
zmax3b = zeros(size(max3));
[i3,j3,~] = ind2sub(size(max3),find(max3));
t3 = cum_max2(max3);
ind3 = sub2ind(size(num_max2),i3,j3,t3);
zmax3b(max3) = num_max2(ind3);

% max3a = max3;
% Ir = 1;
% max3a(:,:,zmin(Ir):zmax(Ir)) = max3a(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
% 
% max3b = max3;
% Ir = 1;
% max3b(:,:,zmin(Ir):zmax(Ir)) = max3b(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
% 
% max3 = max3a | max3b;


Ar3 = -Nr3:Nr3;
zmin3 = max(1,(1+Ar3));
zmax3 = min(N_layer,(N_layer+Ar3));
for Ir = 1:length(Ar3)
    max30(:,:,zmin3(Ir):zmax3(Ir)) = max30(:,:,zmin3(Ir):zmax3(Ir)) | imdilate(max3(:,:,zmin3(end-Ir+1):zmax3(end-Ir+1)),strel('square',Ixy));
end
max30 = max20 & max30;

mask_out = max3;

mask_out = mask_out & (imstack*65535 >= im_max);
max30 = max30 & (imstack*65535 >= im_max);

varargout = {max30 & (~mask_out),zmax3b};