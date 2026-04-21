clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_folder = 'Calibration3_01052020/';
lsm_type = '*.lsm';
out_folder = 'stacks/';
tif_name = 'stack';
figure_tail = '.tif';
match_file = 'matchlist.xls';
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
        lsm_stack = tiffread([input_folder,input_name]); %%% Lsm file loading

        for I_layer = 1:length(lsm_stack)
            tiff_image = lsm_stack(I_layer).data;

        
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
        file_name{I_file,2} = lsm_stack(1).lsm.VoxelSizeX/(1e-6); %%% resolution
%%% =======================================================================
    end

%file_num(strmatch(else_name,{lsm_name.name}),:) = [match_layer,eNbin,eMdim,eWGA_channel,compare_ratio,eWGA_channel,eDAPI_channel,esignal_channel,resolution];
    xlswrite([input_folder,out_folder,match_file],file_name);

