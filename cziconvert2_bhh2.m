function cziconvert2(lf_name)
%clear all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to convert czi file to tiff stacks without tiling %%%%%%%%%%%
%% lf_name (input): txt: excel file names for input data; %%%%%%%%%%%%%%%%%
%%                  cell: input data expressed as cell data. %%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsm_type = '*.czi';
out_folder = 'stacks/';
tile_name = 'tile';
time_name = 'time';
tif_name = 'stack';
figure_tail = '.tif';
match_file = 'matchlist.xls';
compare_ratio = 1;

fast_add = {'Fast','fast'};
airy_add = {'Airy','airy'};
obj60_add = {'60X','60x'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name = lf_name;
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
    pathname=replace(lf_name,'Duallist.xls','');
    folder_new=cellfun(@(x) [pathname,x],folder_list(:,1),'UniformOutput',false);
%     folder_list(:,1)=folder_new;
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
    if size(lsm_name,1)==0
        input_folder = folder_new{list_I,1};
        lsm_name = dir([input_folder,lsm_type]);
    end
    if exist([input_folder,out_folder]) ~= 7
%         mkdir([input_folder,out_folder]);
    end

    file_name = cell(length(lsm_name),3);
    file_num = zeros(length(lsm_name),12+length(signal3_channel)-1);
    output_I = 0;

    for I_file = 1:length(lsm_name)
        input_name = lsm_name(I_file).name;
%         lsm_stack = bfopen([input_folder,input_name]); %%% czi file loading
        reader = bfGetReader([input_folder,input_name]); %%% access czi file
        
        N_C = reader.getSizeC; %%% get the number of colors
        N_Z = reader.getSizeZ; %%% get the number of Z slices
        N_T = reader.getSizeT; %%% get the number of time points
        N_S = reader.getSeriesCount(); %%% get the number of tiles
        N_S0 = prod(Ntile); %%% get the number of tiles per embryo
        N_E = N_S/N_S0; %%% get the number of embryos

        for I_embryo = 1:N_E
            output_I = output_I+1;
%%% Output folder setup: %%% ------------------------------------------
            if  N_E > 1
                pack_name = ['_',num2str(I_embryo,'%03u')];
            else
                pack_name = '';
            end
            
            output_name = [input_name(1:(find(input_name == '.',1,'last')-1)),pack_name,'/'];
            if exist([input_folder,out_folder,output_name]) ~= 7
                mkdir([input_folder,out_folder,output_name]);
            end
%%% -------------------------------------------------------------------

            for I_series = 1:N_S0
                output_name2 = [tile_name,num2str(I_series,'%02u'),'/'];
                if exist([input_folder,out_folder,output_name,output_name2]) ~= 7
                    mkdir([input_folder,out_folder,output_name,output_name2]);
                end
                 reader.setSeries((I_embryo-1)*N_S0+I_series-1);
%                reader.setSeries((I_embryo-1)+(I_series-1)*N_E);

                for I_T = 1:N_T

                    for I_layer = 1:N_Z
                        tiff_image = zeros(0);

                        for I_C = 1:N_C
                            iPlane = reader.getIndex(I_layer-1,I_C-1,I_T-1)+1;
                            image_temp = bfGetPlane(reader,iPlane);
                            tiff_image = cat(3,tiff_image,image_temp);
                        end

                        if size(tiff_image,3) == 1
                            tiff_image(:,:,2) = tiff_image(:,:,1);
                        end
                        if size(tiff_image,3) == 2
                            tiff_image(:,:,3) = tiff_image(:,:,2);
                        end
%                         if size(tiff_image,3) == 3
%                             tiff_image(:,:,4) = tiff_image(:,:,2);
%                         end
%%% Image output: %%% -------------------------------------------------
                        out_num1 = num2str(I_layer,'%02u');
                        out_num2 = num2str(I_T,'%04u');
                        out_stack = [time_name,out_num2,[tif_name,out_num1,figure_tail]];
                        tiffwrite0(tiff_image,[input_folder,out_folder,output_name,output_name2,out_stack]);
%%% -------------------------------------------------------------------
                    end
                end
            end

%%% matchlist.xls output: %%%==============================================
            file_name{output_I,1} = output_name;   %%% tiff folder name
            %%% Data type: F: fast airyscan; A: airyscan; D: double direction tiling
            if ~isempty(strfind(output_name,fast_add{1})) || ~isempty(strfind(output_name,fast_add{2})) || ~isempty(strfind(folder_list{list_I,7},fast_add{1})) || ~isempty(strfind(folder_list{list_I,7},fast_add{2}))
                file_name{output_I,2} = 'FD';
            elseif ~isempty(strfind(output_name,airy_add{1})) || ~isempty(strfind(output_name,airy_add{2})) || ~isempty(strfind(folder_list{list_I,7},airy_add{1})) || ~isempty(strfind(folder_list{list_I,7},airy_add{2}))
                file_name{output_I,2} = 'AD';
            else
                file_name{output_I,2} = 'D';
            end

            if ~isempty(strfind(output_name,obj60_add{1})) || ~isempty(strfind(output_name,obj60_add{2})) || ~isempty(strfind(folder_list{list_I,7},obj60_add{1})) || ~isempty(strfind(folder_list{list_I,7},obj60_add{2}))
                file_name{output_I,2} = [file_name{output_I,2},obj60_add{1}];
            end

            %match_layer = ceil(length(lsm_stack)/2);
            if isempty(match_list)
                match_layer = 1;
            else
                match_layer = match_list(output_I);
            end

            omeMeta = reader.getMetadataStore();
            voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
            resolution = voxelSizeX.doubleValue();
            try
                voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER); % in µm
                resolutionz = voxelSizeZ.doubleValue();        
            catch
                resolutionz = 0;        
            end
            file_num(output_I,:) = [match_layer,Ntile,match_channel,compare_ratio,WGA_channel,DAPI_channel,signal1_channel,resolution,signal2_channel,resolutionz,signal3_channel];
    %         file_num(output_I,:) = [match_layer,Nbin,Mdim,match_channel,compare_ratio,WGA_channel,DAPI_channel,signal1_channel,resolution,signal2_channel,resolutionz,signal3_channel];
%%% =======================================================================
        end

        reader.close()
    end

%file_num(strmatch(else_name,{lsm_name.name}),:) = [match_layer,eNbin,eMdim,eWGA_channel,compare_ratio,eWGA_channel,eDAPI_channel,esignal_channel,resolution];
    try
        xlswrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
    catch
%         xlwrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
        writecell(cat(2,file_name,num2cell(file_num)),[input_folder,out_folder,match_file]);
    end
end
