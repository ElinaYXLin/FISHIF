function delete_mask(threshold_single,threshold_whole,edge,DAPI_channel,AreaRange,maskpath,stackpath)
% threshold=0.3;DAPI_channel=4;
% maskpath='Y:\Imaging\S1\05222020S1_60X\masks\05222020-bcd-488-hb1-TMR-hb7-647_60X_2020_05_22__10_33_53_001_new\';
% stackpath='Y:\Imaging\S1\05222020S1_60X\stacks\05222020-bcd-488-hb1-TMR-hb7-647_60X_2020_05_22__10_33_53_001_new\';
mask_nmae  = dir([maskpath,'*mask.mat']);
time_dim = size(dir([stackpath,'*stack01.tif']),1);
bw_path=[stackpath,'Intermediate\'];
for time_index=1:time_dim
    try
        load([maskpath,'time',num2str(time_index,'%04u'),'mask.mat']);
    end
    stack_name=dir([stackpath,'time',num2str(time_index,'%04u'),'stack*.tif']);
    if size(stack_name,1)==0 && time_dim==1
        stack_name=dir([stackpath,'stack*.tif']);
        load([maskpath,'mask.mat']);
    end
    stack_name={stack_name.name};
    sg_layers=nan(max(mask_stack(:)),size(mask_stack,3));
    AreaDelG={};
    for I_layer = 1:size(mask_stack,3)
        img=imread([stackpath,stack_name{I_layer}]);
        sg_prop = regionprops(double(mask_stack(:,:,I_layer)),double(img(:,:,DAPI_channel)),'MeanIntensity');
        sgArea = regionprops(double(mask_stack(:,:,I_layer)),'Area');
        AreaDelIndex=[sgArea.Area]'>AreaRange(2)|[sgArea.Area]'<AreaRange(1);
        AreaNum=(1:size(sgArea,1))';
        AreaDelG{I_layer}=AreaNum(AreaDelIndex);%del out of area range
        sg_layers(1:size(sg_prop,1),I_layer)=[sg_prop.MeanIntensity]';
    end
    
%     sg_max=max(sg_layers(:));
    sg_label_max=max(sg_layers,[],2);
    sg_label_order=sort(sg_label_max,1);
    sg_middle=sg_label_order(floor(size(sg_label_order,1)/2));
    low_index=sg_label_max<sg_middle*threshold_whole;
    Index_num=1:length(low_index);
    low_index=Index_num(low_index);
    for ii=1:size(low_index,2)
        mask_stack(mask_stack==low_index(ii))=0;
    end
    
    [sg_label_max,max_index]=max(sg_layers,[],2);
    sg_layers_del=sg_layers<sg_label_max*threshold_single;
    sg_layers_notnan=~isnan(sg_layers);
    for ii=1:size(sg_layers,1)
        sg_layers_del_left=sg_layers_del(ii,1:max_index(ii));
        [~,on_index]=find(sg_layers_del_left==1);
        if size(on_index,2)>=1
            sg_layers_del(ii,1:on_index(end))=1;
        end
        sg_layers_del_right=sg_layers_del(ii,max_index(ii):end);
        [~,on_index]=find(sg_layers_del_right==1);
        if size(on_index,2)>=1
            sg_layers_del(ii,(max_index(ii)+on_index(1)-1):end)=1;
        end
    end
    sg_layers_del=sg_layers_del&sg_layers_notnan;
    for I_layer = 1:size(mask_stack,3)
        mask_stack_del=mask_stack(:,:,I_layer);
        %% Area Del
        AreaDel=AreaDelG{I_layer};
        for ii=1:size(AreaDel,1)
            mask_stack_del(mask_stack_del==AreaDel(ii))=0;
        end
        %% intensity Del
%         low_index=sg_layers(:,I_layer)<sg_label_max*threshold_single;
        low_index=sg_layers_del(:,I_layer);
        Index_num=1:length(low_index);
        low_index=Index_num(low_index);
        for ii=1:size(low_index,2)
            mask_stack_del(mask_stack_del==low_index(ii))=0;
        end
        %% edge Del
        if edge==1
%             re_path=strrep(maskpath,'\masks\','\Results\');
            [startIndex,endIndex]=regexp(maskpath,'\masks\');
            re_path=[maskpath(1:startIndex-1),'Results\'];
            re_name=[maskpath(endIndex+1:end-1),'.mat'];
            load([re_path,re_name],'em_mask')
            mask_out=mask_stack_del(~em_mask);
            mask_out=unique(mask_out);
            for ii=1:length(mask_out)
                mask_stack_del(mask_stack_del==mask_out(ii))=0;
            end
        end
        bw_masks=logical(mask_stack_del);
        save([bw_path,'bw_',stack_name{I_layer}(1:end-4),'_cp_masks.mat'],'bw_masks')
    end
end
end
    