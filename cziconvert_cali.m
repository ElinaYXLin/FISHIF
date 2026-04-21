clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input_folder = 'Calibration3_FA_08112019/';
input_folder = 'Calibration_Inten_TMR_FA_01052020b/';
lsm_type = '*.czi';
out_folder = 'stacks/';
tif_name = 'stack';
figure_tail = '.tif';
match_file = 'matchlist.xls';
offset0 = 10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%% LSM file loading/resave: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lsm_name = dir([input_folder,lsm_type]);
    if exist([input_folder,out_folder]) ~= 7
        mkdir([input_folder,out_folder]);
    end

    file_name = cell(length(lsm_name)/2,2);
    output_I = 0;
    output_switch = false;

    for I_file = 1:length(lsm_name)
        input_name = lsm_name(I_file).name;
        lsm_stack = bfopen([input_folder,input_name]); %%% czi file loading

        for I_layer = 1:size(lsm_stack{1,1},1)
            tiff_image = lsm_stack{1,1}{I_layer,1}-offset0;

        
    %%% Image output: %%% -------------------------------------------------
            if length(lsm_stack) > 1
                out_num = num2str(I_layer,'%02u');
            else
                out_num = '';
            end
            out_stack = [input_name(1:(find(input_name == '.',1,'last')-1)),out_num,figure_tail];
%             imwrite(tiff_image,[input_folder,out_folder,out_stack]);
            tiffwrite0(tiff_image,[input_folder,out_folder,out_stack]);
    %%% -------------------------------------------------------------------
        end
    
%%% matchlist.xls output: %%%==============================================
        file_name{I_file,1} = out_stack;
        voxelSizeX = lsm_stack{1,4}.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
        resolution = voxelSizeX.doubleValue();
        file_name{I_file,2} = resolution; %%% resolution
%%% =======================================================================
%         lsm_stack.close()
    end

%file_num(strmatch(else_name,{lsm_name.name}),:) = [match_layer,eNbin,eMdim,eWGA_channel,compare_ratio,eWGA_channel,eDAPI_channel,esignal_channel,resolution];
    xlswrite([input_folder,out_folder,match_file],file_name);

