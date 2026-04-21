function [outimage,IImatch,mismatch] = tilematch(input1,input2,Mdim,match_channel,prematch,compare_ratio)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tiling image matching and grafting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function merges two images that have some overlap on one dimension, and output a bigger image.
%%% outimage: (image name), output merged image.
%%% IImatch: (integer), merge area width.
%%% mismatch: (number 0~1), mismatch value (std).
%%% input1: (image name), input image 1 to merge.
%%% input2: (image name), input image 2 to merge.
%%% Mdim: (integer: 1,2), merge dimension (1: up-down, 2: left-right).
%%% match_channel: (integer: 1,2,3), channel # for matching optimization calculation.
%%% prematch: (integer), pre_setted match displacement. if the value is 0, the program will search for the optimal value, otherwise, the prematch value will be used to merge images.
%%% compare_ratio: (number 0~1), ratio of total pixels (along the non-matching dimension) to use for match optimization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting/intialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mthreshold = 1;
max_width = 400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Image matching and grafting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initiation: %%%========================================================
Nsize1 = size(input1);
Nsize2 = size(input2);
dim_size = length(Nsize1);

if (dim_size > 3) || (dim_size ~= length(Nsize2)) || any(Nsize1(1:dim_size ~= Mdim) ~= Nsize2(1:dim_size ~= Mdim));
    error(['Image size mistake: input1 = ',num2str(Nsize1),'; input2 = ',num2str(Nsize2)])
elseif (dim_size < 3) || (match_channel == 1)
    Nsize1(3) = 1;
    Nsize2(3) = 1;
elseif (dim_size < 3) && (match_channel > 1)
    error(['Image size mistake (single channel image): ',input_folder,input_file1,' = ',num2str(Nsize1),'; ',input_folder,input_file2,' = ',num2str(Nsize2),'; channel = ',num2str(match_channel)])
end

Nmatch = min(Nsize1(Mdim),Nsize2(Mdim)); %%% scanning size
input1_start = uint16([(1-compare_ratio)/2,(1-compare_ratio)/2].*Nsize1(1:2)); %%% scanning start region of input1
input1_end = uint16([(1+compare_ratio)/2,(1+compare_ratio)/2].*Nsize1(1:2)); %%% scanning end region of input1
input2_start = uint16([(1-compare_ratio)/2,(1-compare_ratio)/2].*Nsize2(1:2)); %%% scanning start region of input2
input2_end = uint16([(1+compare_ratio)/2,(1+compare_ratio)/2].*Nsize2(1:2)); %%% scanning end region of input2
input1_start = input1_start+1;
input2_start = input2_start+1;
input1_end(Mdim) = Nsize1(Mdim); %%% size reset for the matching dimension
input2_start(Mdim) = 1; %%% size reset for the matching dimension

mismatch = ones(1,Nmatch); %%% mean mismatch value for all overlaps (original map)
mismatch2 = ones(1,Nmatch); %%% mean mismatch value for all overlaps (bit map)
match_ratio = ones(1,Nmatch); %%% mean intensity ratio
%%%========================================================================

%%% Search the whole region to find the minimal mismatch: %%%==============
if prematch > 0
    IImatch = prematch;
    match_start = IImatch;
    match_end = IImatch;
else
    match_start = 1;
    match_end = Nmatch;
end

for Imatch = match_start:match_end
    input1_start(Mdim) = Nsize1(Mdim)-Imatch+1;
    input2_end(Mdim) = Imatch;
    if Imatch > max_width %%% if the overlap is larger than the max window, put window in the middle of the overlap area
        I_reduce = floor(Imatch/2);
        input1_start(Mdim) = Nsize1(Mdim)-I_reduce-floor(max_width/2)+1;
        input1_end(Mdim) = Nsize1(Mdim)-I_reduce+floor(max_width/2);
        input2_start(Mdim) = I_reduce-floor(max_width/2)+1;
        input2_end(Mdim) = I_reduce+floor(max_width/2);
    end
    temp1 = input1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),match_channel);
    temp2 = input2(input2_start(1):input2_end(1),input2_start(2):input2_end(2),match_channel);
    match_mask = im2bw(uint16(temp1),mthreshold*graythresh(uint16(temp1)))|im2bw(uint16(temp2),mthreshold*graythresh(uint16(temp2)));
    match_ratio(Imatch) = sum(sum(double(temp1)))/sum(sum(double(temp2)));
    raw_mismatch = (double(temp1)-match_ratio(Imatch)*double(temp2))./(double(temp1)+match_ratio(Imatch)*double(temp2)+((double(temp1)+match_ratio(Imatch)*double(temp2))<=0));
    mismatch(Imatch) = sum(sum(raw_mismatch.*raw_mismatch.*double(match_mask)))/sum(sum(match_mask));
    mismatch2(Imatch) = sum(sum(double((im2bw(uint16(temp1),mthreshold*graythresh(uint16(temp1)))-im2bw(uint16(temp2),mthreshold*graythresh(uint16(temp2)))) ~= 0)))/(sum(sum(double(match_mask)))+(sum(sum(double(match_mask)))<=0));
end
[~,IImatch] = min(mismatch);
[~,IImatch2] = min(mismatch2);

if prematch > 0
    IImatch = prematch;
end
%%%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Index set up: %%%======================================================
input1_start = [1,1]; %%% scanning start region of input1 for merge region
input1_end = Nsize1(1:2); %%% scanning end region of input1 for merge region
input2_start = [1,1]; %%% scanning start region of input2 for merge region
input2_end = Nsize2(1:2); %%% scanning end region of input2 for merge region
input1_start(Mdim) = Nsize1(Mdim)-IImatch+1; %%% reset the start point
input2_end(Mdim) = IImatch; %%% reset the end point
if IImatch > max_width
    input1_start(Mdim) = Nsize1(Mdim)-max_width+1;
    input2_start(Mdim) = IImatch-max_width+1;
end

input0_start = [1,1]; %%% scanning start region of input1
input0_end = Nsize1(1:2); %%% scanning end region of input1
input3_start = [1,1]; %%% scanning start region of input2
input3_end = Nsize2(1:2); %%% scanning end region of input2
I_reduce = floor(IImatch/2);
input0_end(Mdim) = Nsize1(Mdim)-I_reduce;
input3_start(Mdim) = I_reduce+1;
%input3_start(Mdim) = IImatch+1; %%% reset the start point
%if IImatch > max_width
%    input0_end(Mdim) = Nsize1(Mdim)-max_width; %%% reset the start point
%else
%    input0_end(Mdim) = Nsize1(Mdim)-IImatch; %%% reset the start point
%end
%%%========================================================================

%%% Image grafting: %%%====================================================
%outimage = uint16(cat(Mdim,input1(input0_start(1):input0_end(1),input0_start(2):input0_end(2),:),(input1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:)+input2(input2_start(1):input2_end(1),input2_start(2):input2_end(2),:))/2,input2(input3_start(1):input3_end(1),input3_start(2):input3_end(2),:)));
%outimage = uint16(cat(Mdim,input1,input2(input3_start(1):input3_end(1),input3_start(2):input3_end(2),:)));
outimage = uint16(cat(Mdim,input1(input0_start(1):input0_end(1),input0_start(2):input0_end(2),:),match_ratio(IImatch)*input2(input3_start(1):input3_end(1),input3_start(2):input3_end(2),:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%