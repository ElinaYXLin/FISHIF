%immask = zeros(size(imstack));
I_layer = 6;
%for I_layer = 1:size(imstack,3)
    
    imtest = imstack(:,:,I_layer);
    immask2 = im2bw(imtest,graythresh(imstack)*0.65);
    immask2 = imclose(immask2,strel('disk',3));
    immask2 = imopen(immask2,strel('disk',20));
    immask2 = bwareaopen(immask2,1000);
    se_mask = regionprops(immask2,'Area');
    se_area = [se_mask.Area];
    bwout = false(size(immask2));   %%% Output mask initialization

    %% Single nucleus area finding: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    area1 = geomean(se_area(se_area > 1000));   %%% calculate the mean area for a single nucleus
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Partition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bwraw0 = immask2;
    bwraw = bwraw0;
    bwraw = bwareaopen(bwraw,round(area1*0.5));   %%% remove small areas
    nthin = 0;
    bthin = 5;

    while any(any(bwraw)) && (nthin <= 50)
        se_mask = regionprops(bwraw,'Area');   %%% area calculation
        nmask = round([se_mask.Area]./area1);   %%% calculate the estimated # of nuclei for each recognized area
        lmask = bwlabel(bwraw);   %%% mask labeling

        bwout(ismember(lmask,find(nmask <= 1))) = true;   %%% Add the single nuclei areas to the output mask
        bwraw(ismember(lmask,find(nmask <= 1))) = false;   %%% Remove the single nuclei areas from the raw mask

        bwraw = imclose(bwraw,strel('disk',3));
        
        nthin = nthin+bthin   %%% thinning value 
        bwthin = bwmorph(bwraw,'skel',nthin);
        %bwthin = imerode(bwraw,strel('disk',nthin));
        %%% raw mask thinning
        bwthin = bwmorph(bwthin,'spur',Inf);   %%% remove spur pixels
        bwthin(ismember(imfilter(uint8(bwthin),ones(3),'symmetric','conv'),[3,2,1])) = false;   %%% remove lines
        bwline = (~bwmorph(bwthin,'thicken',(nthin+100))) & bwraw;
%        bwline = bwconvhull(bwline,'objects');
        tline = bwmorph(bwline,'thicken',1) & bwraw0;

        bwtemp = bwlabel(bwraw & (~bwline));
        se_temp = regionprops(logical(bwtemp),'Area');
        temp_area = [1e9,[se_temp.Area]];
        tempI = temp_area(bwtemp+1);
        se_line = regionprops(tline,tempI,'MinIntensity','MaxIntensity');
        smline = find(([se_line.MinIntensity] < area1*0.5) & ([se_line.MinIntensity] ~= [se_line.MaxIntensity]));
        tline(ismember(bwlabel(tline),smline)) = false;
    %     a(:,:,3) = bwraw;
    %     a(:,:,2) = tline & bwraw;
    %     a(:,:,1) = 0;
    %     imshow(double(a));
    %     title(num2str(nthin))
    %     waitforbuttonpress
        bwraw(tline) = false;   %%% Add new partition boundaries
    end
    %bwout = bwareaopen(bwout,round(area1*0.5));   %%% remove small areas
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    bwout1 = imerode(bwout,strel('disk',3));
    bwout1 = bwconvhull(bwconvhull(bwout1,'objects'),'objects');
    bwout1 = bwmorph(bwout1,'thicken',3);
%    immask(:,:,I_layer) = bwout1;

%end
