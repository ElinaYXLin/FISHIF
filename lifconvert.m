function lifconvert(lf_name)
%clear all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsm_type = '*.lif';
out_folder = 'stacks/';
tif_name = 'stack';
figure_tail = '.tif';
match_file = 'matchlist.xls';
match_key = 'match';
compare_ratio = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name = lf_name;
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);

for list_I = 1:N1
    input_folder = folder_list{list_I,1};
    if isempty(folder_list{list_I,2})
        match_list = [];
    else
        match_list = eval(folder_list{list_I,2});
    end
    all_channel = eval(folder_list{list_I,3});
    DAPI_channel = all_channel(1);
    WGA_channel = all_channel(2);
    signal1_channel = all_channel(3);
    signal2_channel = all_channel(4);
    if length(all_channel) > 4
        signal3_channel = all_channel(5:end);
    else
        signal3_channel = zeros(1);
    end
    match_channel = WGA_channel;
    all_other = eval(folder_list{list_I,4});
    Nbin = all_other(1);
    Mdim = all_other(2);
    
%% LSM file loading/resave: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lsm_name = dir([input_folder,lsm_type]);
    if exist([input_folder,out_folder]) ~= 7
        mkdir([input_folder,out_folder]);
    end

    file_name = cell(length(lsm_name),2);
    file_num = zeros(length(lsm_name),12+length(signal3_channel)-1);
    output_I = 0;
    output_switch = false;

    for I_file = 1:length(lsm_name)
        input_name = lsm_name(I_file).name;
%         lsm_stack = tiffread([input_folder,input_name]); %%% Lsm file loading
        [imgout]=ci_loadLif([input_folder,input_name]); %%% Lif file loading
        lsm_stack = lif_to_stack(imgout);
        
        pack_num = 1;
        stack_num = 0;
        bin_num = 0;
        if  ~isempty(strfind(input_name,match_key))
            bin_size = Nbin*2-1+(Nbin == 1);
        else
            bin_size = 1;
        end
        if ~(ismember(input_name(find(input_name == '.',1,'last')-1),'Bb') && input_name(find(input_name == '.',1,'last')-2) == '_')
            output_I = output_I+1;
        end
        
        while stack_num(end) < length(lsm_stack)
            bin_num = bin_num+1;
            pack_num = pack_num+(bin_num > bin_size);
            output_I = output_I+(bin_num > bin_size);
            bin_num = mod(bin_num,bin_size);
            if bin_num == 0
                bin_num = bin_size;
            end
            stack_num = stack_num(end)+[1:lsm_stack(stack_num(end)+1).lsm.DimensionZ];
            
    %%% Output folder setup: %%% ------------------------------------------
            if  bin_size*lsm_stack(1).lsm.DimensionZ >= length(lsm_stack)
                pack_name = '';
            else
                pack_name = ['_',num2str(pack_num,'%03u')];
            end
            if ~isempty(strfind(input_name,match_key))
                match_name = ['_',char('B'-mod(bin_num,2))];
            else
                match_name = '';
            end
            output_name = [input_name(1:(find(input_name == '.',1,'last')-1)),pack_name,match_name,'/'];
            if exist([input_folder,out_folder,output_name]) ~= 7
                mkdir([input_folder,out_folder,output_name]);
            end
    %%% -------------------------------------------------------------------
    
            for I_layer = stack_num
                stack_raw = lsm_stack(I_layer).data;
                tiff_image = zeros(0);
                if iscell(stack_raw)
                    for I_color = 1:length(stack_raw)
                        tiff_image = cat(3,tiff_image,stack_raw{I_color});
                    end
                else
                    tiff_image = stack_raw;
                end
                
                if size(tiff_image,3) == 1
                    tiff_image(:,:,2) = tiff_image(:,:,1);
                end
                if size(tiff_image,3) == 2
                    tiff_image(:,:,3) = tiff_image(:,:,2);
                end
    %%% Image output: %%% -------------------------------------------------
                if length(lsm_stack) > 1
                    out_num = num2str(I_layer-stack_num(1)+1,'%02u');
                else
                    out_num = '';
                end
                out_stack = [tif_name,out_num,figure_tail];
                if bin_num <= 2
%                     imwrite(tiff_image,[input_folder,out_folder,output_name,out_stack]);
                    tiffwrite0(tiff_image,[input_folder,out_folder,output_name,out_stack]);
                else
                    temp_old = imread([input_folder,out_folder,output_name,out_stack]);
                    tiff_image = cat(Mdim,temp_old,tiff_image);
%                     imwrite(tiff_image,[input_folder,out_folder,output_name,out_stack]);
                    tiffwrite0(tiff_image,[input_folder,out_folder,output_name,out_stack]);
                end
    %%% -------------------------------------------------------------------
            end
    
%%% matchlist.xls output: %%%==============================================
            if ismember(output_name(length(output_name)-1),'Bb') && output_name(length(output_name)-2) == '_'
                output_switch = true;
            else
                output_switch = false;
            end
            file_name{output_I,output_switch+1} = output_name;
            %match_layer = ceil(length(lsm_stack)/2);
            if isempty(match_list)
                match_layer = 1;
            else
                match_layer = match_list(output_I);
            end
            
            resolution = lsm_stack(1).lsm.VoxelSizeX/(1e-6);
            resolutionz = lsm_stack(1).lsm.VoxelSizeZ/(1e-6);
            file_num(output_I,:) = [match_layer,Nbin,Mdim,match_channel,compare_ratio,WGA_channel,DAPI_channel,signal1_channel,resolution,signal2_channel,resolutionz,signal3_channel];
%%% =======================================================================
        end
    end

%file_num(strmatch(else_name,{lsm_name.name}),:) = [match_layer,eNbin,eMdim,eWGA_channel,compare_ratio,eWGA_channel,eDAPI_channel,esignal_channel,resolution];
    try
        xlswrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
    catch
        xlwrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
    end
end



function lsm_stack = lif_to_stack(imgout)

N_layer = size(imgout.Image{1},3);
N_color = size(imgout.Image,1);

for I_layer = 1:N_layer
    lsm_stack(I_layer).lsm.DimensionZ = N_layer;
    
    lsm_stack(I_layer).data = cell(N_color,1);
    for I_color = 1:N_color
        lsm_stack(I_layer).data{I_color} = imgout.Image{I_color}(:,:,I_layer);
    end
    
    [Lx,Ly,Lz] = imgout.Info.Dimensions.Length;
    [Nx,Ny,Nz] = imgout.Info.Dimensions.NumberOfElements;
    lsm_stack(I_layer).lsm.VoxelSizeX = abs(str2double(Lx)/str2double(Nx));
    lsm_stack(I_layer).lsm.VoxelSizeZ = abs(str2double(Lz)/str2double(Nz));

end



