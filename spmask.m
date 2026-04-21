function [mask_out,varargout] = spmask(imstack,Nr,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% A function for spot recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mask_out: Output spot recognition result;
%% imstack: input image stack;
%% Nr: spot spread half width (on z direction);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max3 = imregionalmax(imstack);   %%% 3D local maxima
max2 = false(size(max3));   %%% 2D local maxima
max20 = max2;
N_layer = size(imstack,3);
for I_layer = 1:N_layer
    max20(:,:,I_layer) = imregionalmax(imstack(:,:,I_layer));
    max2(:,:,I_layer) = imdilate(max20(:,:,I_layer),strel('square',3));
end

Nr = ceil(Nr);
Ar = -Nr:Nr;
zmin = max(1,(1+Ar));
zmax = min(N_layer,(N_layer+Ar));
for Ir = 1:length(Ar)
    max3(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
end

mask_out = max3;

if ~isempty(varargin)
    mask_out = mask_out & (imstack*65535 >= varargin{1});
    max20 = max20 & (imstack*65535 >= varargin{1});
end

varargout = {max20 & (~mask_out)};