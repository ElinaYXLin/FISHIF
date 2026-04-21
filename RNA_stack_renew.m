function RNA_stack_renew(imfolder_list,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to refit foci %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% imfolder_list: list of data files.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_name = 'matchlist.xls';
if isempty(varargin)
    input_mark = 'T';
else
    input_mark = varargin{1};
end

if isa(imfolder_list,'char') && strcmp(imfolder_list(end-3:end),'.xls')
    list_name0 = imfolder_list;
    [~, txt_list] = xlsread(list_name0);
    imfolder_list = txt_list(strcmpi(input_mark,txt_list(:,6)),1);
elseif isa(imfolder_list,'char')
    error('Incorrect input!')
end

nd_all = cell(0);
for I_folder = 1:length(imfolder_list)
    imfolder = imfolder_list{I_folder};
    if ~(imfolder(end) == '/' || imfolder(end) == '\')
        imfolder = [imfolder,'/'];
    end
    
    I_cut = strfind(imfolder,'stacks');
    if I_cut
        imname = imfolder(I_cut+7:end-1);
        imfolder = imfolder(1:I_cut+6);
        n0 = foci_renew(imfolder,imname);
        nd_all = cat(1,nd_all,{imfolder,imname,n0});
    else           
        imfolder = [imfolder,'stacks/'];
        [~,file_list] = xlsread([imfolder,list_name]);
        if size(file_list,1) > 0
            for I_file = 1:size(file_list,1)
                imname = file_list{I_file,3};
                imname = imname(1:(end-1));
                n0 = foci_renew(imfolder,imname);
                nd_all = cat(1,nd_all,{imfolder,imname,n0});
            end
        end
    end
    
end

nd_all
xlswrite('RNA_stack_renew_log.xls',nd_all)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function n_diff = foci_renew(imfolder,imname)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to refit foci for a particular data file %%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% n_diff (output): foci number difference from the previous result;
%% imfolder (input): image set folder;
%% imname (input): image name.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'matchlist.xls';
image_tail = '*.tif';
hist_folder0 = 'Histogram/';
hist_add = '_raw';
out_tail = '.xls';
out2_add = '_value';
out2_tail = '.mat';
rxy = 10;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image stack loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, file_list] = xlsread([imfolder,list_name]);
I_file = strcmp([imname,'/'],file_list(:,3));
I_file = I_file | strcmp([imname,'/'],file_list(:,3));

lsm_stack = dir([imfolder,imname,'/',image_tail]); %%% image file list loading
RNA_channel = num_list(I_file,10);
immax = length(lsm_stack);
imstack0 = zeros(0);

for I_layer = 1:immax
    temp0 = imread([imfolder,imname,'/',lsm_stack(I_layer).name]);
    imstack0(:,:,I_layer) = temp0(:,:,RNA_channel);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Existing foci list loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_cut = strfind(imfolder,'stacks');
hist_folder = [imfolder(1:I_cut-1),hist_folder0];
hist_name = [imname,hist_add,out_tail];

foci_list = xlsread([hist_folder,hist_name]);
foci_xy = max(round(foci_list(:,6:8)),1);
id_list = sub2ind(size(imstack0),foci_xy(:,1),foci_xy(:,2),foci_xy(:,3));
mask_out = false(size(imstack0));
mask_out(id_list) = true;

foci_2D = any(mask_out,3);
mask_out2D = repmat(foci_2D,[1,1,immax]);
mask_out2D = mask_out2D & (~mask_out);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Foci refitting and output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peak_parameter = ptrack(imstack0,imstack0,mask_out,mask_out2D);   %%% Track/fit mRNA spots
max_spot = spfilter(peak_parameter,rxy,rxy,1);   %%% Fine filter to get rid of noise

xlswrite([hist_folder,hist_name],max_spot);

out2_name = [imname,out2_add,out2_tail];
if exist([hist_folder,out2_name],'file')
    save([hist_folder,out2_name],'max_spot','-append')
end

n_diff = size(max_spot,1)-size(foci_list,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



