function DelPieces(mask_path,varargin)
if ~isempty(varargin) && ~isempty(varargin{1})
    timeIndex = varargin{1};
    mask_name_list=dir([mask_path,'time',num2str(timeIndex,'%04u'),'*mask.mat']);
else
    mask_name_list=dir([mask_path,'*mask.mat']);
end
% mask_name_list=dir([mask_path,'*mask.mat']);
% if size(mask_name_list,1)<=2
%     live_em=0;
% else
%     live_em=1;
% end

bw_path=[strrep(mask_path,'\masks\','\stacks\'),'Intermediate\'];
startI=1;
if strcmp(mask_name_list(1).name,'mask.mat')
    startI=2;
end
for mask_num=startI:size(mask_name_list,1)
    mask_name=mask_name_list(mask_num).name;
    load([mask_path,mask_name]);
%     mask_new=zeros(size(mask_stack),'uint16');
    mkdir([mask_path,'Orginial\'])
    save([mask_path,'Orginial\',mask_name],'mask_stack','-v7.3')
    for I_layer = 1:size(mask_stack,3)
        mask_layer=mask_stack(:,:,I_layer);
        mask_edge=[mask_layer([1 end],:);mask_layer(:,[1 end])'];
        mask_label=unique(mask_edge(:));
        Isedge=~ismember(mask_layer,mask_label);
        mask_new=uint16(double(mask_layer).*Isedge);
        bw_masks=logical(mask_new);
        save([bw_path,'bw_',strrep(mask_name,'mask.mat','stack'),num2str(I_layer,'%02u'),'_cp_masks.mat'],'bw_masks')
    end
end
cellpose_bwmat_live_while(bw_path,mask_path,3)
end