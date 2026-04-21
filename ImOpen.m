function ImOpen(mask_path,DiskSize,Turn,Model)
mask_name_list=dir([mask_path,'*mask.mat']);
% if size(mask_name_list,1)<=2
%     live_em=0;
% else
%     live_em=1;
% end
% Model='ED';
for mask_num=1:size(mask_name_list,1)
    mask_name=mask_name_list(mask_num).name;
    load([mask_path,mask_name]);
    mask_imopen=zeros(size(mask_stack),'uint16');
    mkdir([mask_path,'Orginial\'])
    save([mask_path,'Orginial\',mask_name],'mask_stack','-v7.3')
    for I_layer = 1:size(mask_stack,3)
        mask_layer=mask_stack(:,:,I_layer);
        se = strel('disk',DiskSize);
%         Turn=5;
        afterOpening=mask_layer;
        
        if Model==1
            % imerode-imdilate
            Repet=0;
            while Repet<Turn
                Repet=Repet+1;
                afterOpening = imerode(afterOpening,se);
            end

            Repet=0;
            while Repet<Turn
                Repet=Repet+1;
                afterOpening = imdilate(afterOpening,se);
            end
        else
            Repet=0;
            while Repet<Turn
                Repet=Repet+1;
                afterOpening = imdilate(afterOpening,se);
            end

            Repet=0;
            while Repet<Turn
                Repet=Repet+1;
                afterOpening = imerode(afterOpening,se);
            end
        end
        
%         Repet=0;
%         while Repet<Turn
%             Repet=Repet+1;
%             afterOpening = imopen(afterOpening,se);
%         end
%         figure
%         subplot(1,2,1)
%         imshow(mask_layer,[]);
%         subplot(1,2,2)
%         imshow(afterOpening,[]);
        mask_imopen(:,:,I_layer)=afterOpening;
    end
    mask_stack=mask_imopen;
    save([mask_path,mask_name],'mask_stack','-v7.3')
end
    
end