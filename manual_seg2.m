clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNAlist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
output_tail = '.xls';
figure_tail = '.fig';
mat_tail = '.mat';
seg_add = '_seg';
pro_add = '_pro';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
ex_win = 200;
resolution0 = 0.09;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Manual nuclei segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    for list_J = 1:M1
        
%%% Loading of automatic segmentation result and : %%%==========================
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        L_ratio = (resolution/resolution0);
        imsize = size(seg_bw);
%%% =======================================================================
        
%%% Region selection: %%%==================================================
            x1w = 1;
            x2w = imsize(1);
            y1w = 1;
            y2w = imsize(2);
            local_WGA = max_image(x1w:x2w,y1w:y2w,WGA_channel);
            local_seg = seg_bw(x1w:x2w,y1w:y2w);
            figure(2)
                overlay = imoverlay(adapthisteq(local_WGA), local_seg, [0,0.5,0]);
                imshow(overlay)   %%% show the embryo boundary recognition
                title([image_folder,': segmentation result overview, please select the operation (to add regions, press "a"; to delete regions, press "d"; press "b" when finish)']);
            [~,~,oper_value] = ginput(1);
            while all(oper_value ~= 'Bb')
%%% =======================================================================

%%% Operation: adding regions: %%%=========================================
                if any(oper_value == 'Aa')
                    add_switch = true;
                    while add_switch
                        figure(2)
                            overlay = imoverlay(adapthisteq(local_WGA), local_seg, [0,0.5,0]);
                            title([image_folder,': please click mouse to redraw boundaries, press "del" to go back']);
                        bw_add = roipoly(overlay);
                        add_switch = any(any(bw_add));
                        if add_switch
                            local_seg = local_seg | bw_add;
                        end
                    end
%%% =======================================================================

%%% Operation: adding regions: %%%=========================================
                elseif any(oper_value == 'Dd')
                    add_switch = true;
                    while add_switch
                        figure(2)
                            overlay = imoverlay(adapthisteq(local_WGA), local_seg, [0.5,0,0]);
                            title([image_folder,': please click mouse to redraw boundaries, press "del" to go back']);
                        bw_add = roipoly(overlay);
                        add_switch = any(any(bw_add));
                        if add_switch
                            label_seg = bwlabel(local_seg);
                            local_seg(label_seg == label_seg(find(bw_add,1))) = 0;
                        end
                    end

%%% =======================================================================
                end
                figure(2)
                    overlay = imoverlay(adapthisteq(local_WGA), local_seg, [0,0.5,0]);
                    imshow(overlay)   %%% show the embryo boundary recognition
                    title([image_folder,': Region of Interest, please select the operation (to add regions, press "a"; to delete regions, press "d"; press "b" to go back to the overview)']);
                [~,~,oper_value] = ginput(1);
            end
            
%%% Modification storage: %%%==============================================
            seg_bw(x1w:x2w,y1w:y2w) = local_seg;

%%% =======================================================================

%%% Cytoplasm region recognition: %%%======================================
        embryo_region = false(size(seg_bw));
        if any(any(seg_bw))
            temp_label = 2*seg_bw;
            temp_label(1,1) = 1;
            embryo_prop = regionprops(temp_label,'ConvexImage','BoundingBox');
            x1 = uint16(embryo_prop(2).BoundingBox(2));
            x2 = uint16(x1+embryo_prop(2).BoundingBox(4)-1);
            y1 = uint16(embryo_prop(2).BoundingBox(1));
            y2 = uint16(y1+embryo_prop(2).BoundingBox(3)-1);
            embryo_region(x1:x2,y1:y2) = embryo_prop(2).ConvexImage;
        end
        cyto_bw = embryo_region & (~seg_bw);
        cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1)
        bw_perim_WGA = bwperim(seg_bw);
        %bw_perim_DAPI = bwperim(cumulated_DAPI);
        overlay = imoverlay(adapthisteq(max_image(:,:,WGA_channel)), bw_perim_WGA, [0,1,0]);
        %overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
        imshow(overlay)   %%% show the embryo boundary recognition
        title([image_folder,' (green: nucleus recognition)']);
        %waitforbuttonpress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Result output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),seg_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'seg_bw','cyto_bw','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
end


