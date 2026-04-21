function [signal_stack,RNA_stack,varargout] = stack3D(seg_bw,signal_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch,varargin)

image_type = '*.tif';
imlist = dir([image_folder,image_type]); %%% get the image list from image folder
N_layer = length(imlist);
DAPI_stack = zeros([size(seg_bw),N_layer],'uint16');
signal_stack = zeros([size(seg_bw),N_layer],'uint16');
RNA_stack = zeros([size(seg_bw),N_layer],'uint16');

if ~isempty(varargin)
    signal2_channel = varargin{1};
    protein_signal2_mismatch = varargin{2};
    signal2_stack = zeros([size(seg_bw),N_layer]);
end


%%% Image stack loading: %=================================================
for I_layer = 1:N_layer
    temp = imread([image_folder,imlist(I_layer).name]);
    DAPI_stack(:,:,I_layer) = temp(:,:,DAPI_channel);
    signal_stack(:,:,I_layer) = temp(:,:,signal_channel);
    IR_layer = I_layer+protein_RNA_mismatch;
    if IR_layer <= N_layer && IR_layer > 0
        RNA_stack(:,:,IR_layer) = temp(:,:,RNA_channel);
    end
    if ~isempty(varargin)
    IR_layer2 = I_layer+protein_signal2_mismatch;
        if IR_layer2 <= N_layer && IR_layer2 > 0
            signal2_stack(:,:,IR_layer2) = temp(:,:,signal2_channel);
        end
    end
end

if protein_RNA_mismatch > 0
    RNA_stack(:,:,[1:protein_RNA_mismatch]) = repmat(RNA_stack(:,:,1),[1,1,protein_RNA_mismatch]);
elseif protein_RNA_mismatch < 0
    RNA_stack(:,:,[(end-protein_RNA_mismatch+1):end]) = repmat(RNA_stack(:,:,end),[1,1,protein_RNA_mismatch]);
end

if ~isempty(varargin)
    if protein_signal2_mismatch > 0
        signal2_stack(:,:,[1:protein_signal2_mismatch]) = repmat(signal2_stack(:,:,1),[1,1,protein_signal2_mismatch]);
    elseif protein_signal2_mismatch < 0
        signal2_stack(:,:,[(end-protein_signal2_mismatch+1):end]) = repmat(signal2_stack(:,:,end),[1,1,protein_signal2_mismatch]);
    end
end
%%% =======================================================================

%%% 3D mask generation: ===================================================

%%% =======================================================================
if ~isempty(varargin)
    varargout = {signal2_stack,DAPI_stack};
else
    varargout = {DAPI_stack};
end