%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fly embyro image tiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%folder_name = '01222011/stacks/';
list_name0 = 'GFPlist.xls';
list_sub = 'stacks/';
list_name = 'matchlist.xls';
standard_record = 'Calibration/Results/standard.mat';
output_tail = '.xls';
image_tail = '*.tif';
figure_tail = '.fig';
th_parameter1 = 1.2;
th_parameter2 = 1.2;

global standard_data
standard_data = load(standard_record);

[num_list0, folder_list] = xlsread(list_name0);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
for list_I0 = 1:N1
    folder_name = [folder_list{list_I0,1},list_sub];
    channel_name = eval(folder_list{list_I0,5});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [num_list, file_list] = xlsread([folder_name,list_name]);
    [N1,N2] = size(file_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for list_I = 1:N1
        resolution = num_list(list_I,9);
        if ~isempty(file_list{list_I,2})
    %%% Image list loading: %%%================================================
            imlist1 = dir([folder_name,file_list{list_I,1},image_tail]); %%% get the image list from image folder1
            imlist2 = dir([folder_name,file_list{list_I,2},image_tail]); %%% get the image list from image folder2
            out_name = [file_list{list_I,2}(1:(length(file_list{list_I,2})-3)),'/']; %%% output folder name
            if length(imlist1) ~= length(imlist2)
                error(['Incompatible image stacks: ',folder_name,file_list{list_I,1},' and ',folder_name,file_list{list_I,2}])
            elseif length(imlist1) < num_list(list_I,1)
                error(['Not enough image layers (',folder_name,file_list{list_I,1},'): ',num2str(num_list(list_I,1)),'/',num2str(length(imlist1))])
            end
        %%% =======================================================================

    %%% Initialization of image test matching process: %%%=====================
            Mdim = num_list(list_I,3);
            match_channel = num_list(list_I,4);
            compare_ratio = num_list(list_I,5);
            prematch = 0;
    %%% =======================================================================

    %%% Image loading: %%%=====================================================
            image1 = imread([folder_name,file_list{list_I,1},imlist1(num_list(list_I,1)).name]);
            image2 = imread([folder_name,file_list{list_I,2},imlist2(num_list(list_I,1)).name]);
            Nbin1 = num_list(list_I,2);
            pixel_bin1 = size(image1,Mdim)/Nbin1;
            Nbin2 = Nbin1-1+(Nbin1 == 1);
            pixel_bin2 = size(image2,Mdim)/Nbin2;
        %%%%% Index setup/initialization: %%%%%--------------------------------
            IImatch = zeros(1,(Nbin1+Nbin2-1));
            input1_start = [1,1];
            input1_end = size(image1);
            input1_end(Mdim) = pixel_bin1;
            input2_start = [1,1];
            input2_end = size(image2);
            input2_end(Mdim) = pixel_bin2;
            outimage = image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:);
        %%%%% -----------------------------------------------------------------

        %%%%% Matching parameter search: %%%%%---------------------------------
            Ibin1 = 1;
            Ibin2 = 1;
            I_tile = 1;
            while (Ibin1 < Nbin1)||(Ibin2 <= Nbin2)
                if Ibin2 <= Nbin2
                    input2 = image2(input2_start(1):input2_end(1),input2_start(2):input2_end(2),:);
                    [outimage,IImatch(I_tile),mismatch] = tilematch(outimage,input2,Mdim,match_channel,prematch,compare_ratio);
                    I_tile = I_tile+1;
                    Ibin2 = Ibin2+1;
                    input2_start(Mdim) = input2_start(Mdim)+pixel_bin2;
                    input2_end(Mdim) = input2_end(Mdim)+pixel_bin2;
                end

                if Ibin1 < Nbin1
                    input1_start(Mdim) = input1_start(Mdim)+pixel_bin1;
                    input1_end(Mdim) = input1_end(Mdim)+pixel_bin1;
                    input2 = image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:);
                    [outimage,IImatch(I_tile),mismatch] = tilematch(outimage,input2,Mdim,match_channel,prematch,compare_ratio);
                    I_tile = I_tile+1;
                    Ibin1 = Ibin1+1;
                end
            end
            %IImatch
        %%%%% -----------------------------------------------------------------
    %%% =======================================================================

    %%% Tiling loops: %%%%%====================================================
            for I_layer = 1:length(imlist1)
                image1 = imread([folder_name,file_list{list_I,1},imlist1(I_layer).name]);
                image2 = imread([folder_name,file_list{list_I,2},imlist2(I_layer).name]);
                %%%%% Initialization:
                input1_start = [1,1];
                input1_end = size(image1);
                input1_end(Mdim) = pixel_bin1;
                input2_start = [1,1];
                input2_end = size(image2);
                input2_end(Mdim) = pixel_bin2;
                outimage = corr_dist(image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:),channel_name,resolution);

                %%%%% Tiling:
                Ibin1 = 1;
                Ibin2 = 1;
                I_tile = 1;
                while (Ibin1 < Nbin1)||(Ibin2 <= Nbin2)
                    if Ibin2 <= Nbin2
                        input2 = corr_dist(image2(input2_start(1):input2_end(1),input2_start(2):input2_end(2),:),channel_name,resolution);
                        [outimage,IImatch(I_tile),mismatch] = tilematch2(outimage,input2,Mdim,match_channel,IImatch(I_tile),compare_ratio);
                        I_tile = I_tile+1;
                        Ibin2 = Ibin2+1;
                        input2_start(Mdim) = input2_start(Mdim)+pixel_bin2;
                        input2_end(Mdim) = input2_end(Mdim)+pixel_bin2;
                    end

                    if Ibin1 < Nbin1
                        input1_start(Mdim) = input1_start(Mdim)+pixel_bin1;
                        input1_end(Mdim) = input1_end(Mdim)+pixel_bin1;
                        input2 = corr_dist(image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:),channel_name,resolution);
                        [outimage,IImatch(I_tile),mismatch] = tilematch2(outimage,input2,Mdim,match_channel,IImatch(I_tile),compare_ratio);
                        I_tile = I_tile+1;
                        Ibin1 = Ibin1+1;
                    end
                end

                %%%%% Image output:
                %imshow(outimage);
                %title([folder_name,out_name,': ',num2str(I_layer),'/',num2str(length(imlist1))])
                %waitforbuttonpress
                if exist([folder_name,out_name]) ~= 7
                    mkdir([folder_name,out_name]);
                end
%                 imwrite(outimage,[folder_name,out_name,imlist1(I_layer).name]);
                tiffwrite0(outimage,[folder_name,out_name,imlist1(I_layer).name]);
            end
    %%% =======================================================================

        else
            out_name = file_list{list_I,1}; %%% output folder name
        end

        file_list(list_I,3) = {out_name}; %%% Output folder name record
    end

    xlswrite([folder_name,list_name],cat(2,file_list,num2cell(num_list)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

clear global