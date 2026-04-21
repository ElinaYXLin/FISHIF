function cziconvert(lf_name)
%clear all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsm_type = '*.czi';
out_folder = 'stacks/';
tif_name = 'stack';
figure_tail = '.tif';
match_file = 'matchlist.xls';
match_key = 'match';
overlay_key = 'overlap';
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
    Nbin = all_other(:,1)';
    Mdim = all_other(:,2)';
    
    Ntile = ones(1,2);
    Ntile(Mdim) = Nbin;
    
%% LSM file loading/resave: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lsm_name = dir([input_folder,lsm_type]);
    if exist([input_folder,out_folder]) ~= 7
        mkdir([input_folder,out_folder]);
    end

    file_name = cell(length(lsm_name),3);
    file_num = zeros(length(lsm_name),12+length(signal3_channel)-1);

    for I_file = 1:length(lsm_name)
        input_name = lsm_name(I_file).name;
        lsm_stack = bfopen([input_folder,input_name]); %%% czi file loading
        
        N_C = lsm_stack{1,4}.getChannelCount(0); %%% get the number of colors
        N_Z = lsm_stack{1,4}.getPixelsSizeZ(0).getValue(); %%% get the number of colors
        N_X = Ntile(2);
        N_Y = Ntile(1);

%%% Output folder setup: %%% ------------------------------------------
        output_name = [input_name(1:(find(input_name == '.',1,'last')-1)),'/'];
        if exist([input_folder,out_folder,output_name]) ~= 7
            mkdir([input_folder,out_folder,output_name]);
        end
%%% -------------------------------------------------------------------

        for I_layer = 1:N_Z
            tiff_image = zeros(0);
            
            for I_Y = 1:N_Y
                tiff_temp = zeros(0);
                
                for I_X = 1:N_X
                    tiff0 = zeros(0);
                    
                    for I_C = 1:N_C
                        tiff0 = cat(3,tiff0,lsm_stack{sub2ind([N_X,N_Y],I_X,I_Y),1}{(I_layer-1)*N_C+I_C,1});
                    end
                    
                    tiff_temp = cat(2,tiff_temp,tiff0);
                end
                
                tiff_image = cat(1,tiff_image,tiff_temp);
            end
            
            if size(tiff_image,3) == 1
                tiff_image(:,:,2) = tiff_image(:,:,1);
            end
            if size(tiff_image,3) == 2
                tiff_image(:,:,3) = tiff_image(:,:,2);
            end
%%% Image output: %%% -------------------------------------------------
            out_num = num2str(I_layer,'%02u');
            out_stack = [tif_name,out_num,figure_tail];
            tiffwrite0(tiff_image,[input_folder,out_folder,output_name,out_stack]);
%%% -------------------------------------------------------------------
        end

%%% matchlist.xls output: %%%==============================================
        file_name{I_file,1} = output_name;
        file_name{I_file,1} = 'D';
        %match_layer = ceil(length(lsm_stack)/2);
        if isempty(match_list)
            match_layer = 1;
        else
            match_layer = match_list(I_file);
        end
        
        voxelSizeX = lsm_stack{1,4}.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
        resolution = voxelSizeX.doubleValue();
        voxelSizeZ = lsm_stack{1,4}.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER); % in µm
        resolutionz = voxelSizeZ.doubleValue();        
        
        file_num(I_file,:) = [match_layer,Nbin,match_channel,compare_ratio,WGA_channel,DAPI_channel,signal1_channel,resolution,signal2_channel,resolutionz,signal3_channel];
%         file_num(output_I,:) = [match_layer,Nbin,Mdim,match_channel,compare_ratio,WGA_channel,DAPI_channel,signal1_channel,resolution,signal2_channel,resolutionz,signal3_channel];
%%% =======================================================================
%         lsm_stack.close()
    end

%file_num(strmatch(else_name,{lsm_name.name}),:) = [match_layer,eNbin,eMdim,eWGA_channel,compare_ratio,eWGA_channel,eDAPI_channel,esignal_channel,resolution];
    try
        xlswrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
    catch
        xlwrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
    end
end
