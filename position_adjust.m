function xyz2 = position_adjust(xyz1,max_image,protein_RNA_mismatch,protein_RNA_xymismatch,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to adjust the postion of embryo foci against chromatic aberration: %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% xyz2 (output): adjusted foci positions;
%% xyz1 (input): original foci positions;
%% max_image (input): maximal projection image matrix;
%% protein_RNA_mismatch (input): z layer mismatch between protein and RNA channels
%% protein_RNA_xymismatch (input): xy mismatch between protein and RNA channels (M_protein_RNA,C_protein_RNA,x0_protein_RNA,Nbin,Mdim,resolution,resolution_mismatch)
%% varargin (input): {x_start_im,y_start_im,x_center_im,y_center_im}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xyz2 = zeros(size(xyz1));
xyz2(:,3) = xyz1(:,3)+protein_RNA_mismatch;

%%% chromatic aberration compensation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(protein_RNA_xymismatch)
    M_protein_RNA = protein_RNA_xymismatch{1};
    C_protein_RNA = protein_RNA_xymismatch{2};
    x0_protein_RNA = protein_RNA_xymismatch{3};
    Nbin = protein_RNA_xymismatch{4};
    Mdim = protein_RNA_xymismatch{5};
    resolution = protein_RNA_xymismatch{6};
    resolution_mismatch = protein_RNA_xymismatch{7};
    XY_raw = xyz1(:,1:2);
    
    if  isempty(varargin)
        dim_div = size(max_image,Mdim)/Nbin;
        window_size = [size(max_image,1),size(max_image,2)];
        window_size(Mdim) = dim_div;

        XY_int = zeros(size(XY_raw));
        XY_int(:,Mdim) = floor(XY_raw(:,Mdim)/dim_div)*dim_div;
        XY_mod = XY_raw;
        XY_mod(:,Mdim) = mod(XY_raw(:,Mdim),dim_div);
        XY_center = repmat((1+window_size)/2,size(XY_mod,1),1);
        XY_mod_new = (M_protein_RNA*(XY_mod-XY_center)'+repmat(C_protein_RNA,1,size(XY_mod,1))*resolution_mismatch/resolution)'+XY_center;
        xyz2(:,1:2) = XY_mod_new+XY_int;
    else
        x_start_im = varargin{1}{1}; x_start = zeros(size(x_start_im)); x_end = x_start;
        y_start_im = varargin{1}{2}; y_start = zeros(size(y_start_im)); y_end = y_start;
        x_center_im = varargin{1}{3};
        y_center_im = varargin{1}{4};
        
        x_start = x_start_im-0.5; x_start(:,1) = -inf; x_start3D = repmat(x_start,[1,1,size(XY_raw,1)]);
        x_end(:,1:end-1) = x_start(:,2:end); x_end(:,end) = inf; x_end3D = repmat(x_end,[1,1,size(XY_raw,1)]);
        y_start = y_start_im-0.5; y_start(1,:) = -inf; y_start3D = repmat(y_start,[1,1,size(XY_raw,1)]);
        y_end(1:end-1,:) = y_start(2:end,:); y_end(end,:) = inf; y_end3D = repmat(y_end,[1,1,size(XY_raw,1)]);
        
        y00 = repmat(reshape(XY_raw(:,1),1,1,size(XY_raw,1)),[size(x_start),1]);
        x00 = repmat(reshape(XY_raw(:,2),1,1,size(XY_raw,1)),[size(y_start),1]);
        
        id0 = find(x00 > x_start3D & x00 <= x_end3D & y00 > y_start3D & y00 <= y_end3D);
        [i00,j00,k00] = ind2sub(size(x00),id0);
        [k00,idx] = sort(k00);i00 = i00(idx); j00 = j00(idx);
        id1 = sub2ind(size(x_start),i00,j00);
        if ~isempty(setdiff(k00,1:size(XY_raw,1)))
            error('Fail to adjust spot position: spot does not belong to any tile')
        end
        XY_center = [reshape(y_center_im(id1),size(id1)),reshape(x_center_im(id1),size(id1))];
        
        xyz2(:,1:2) = (XY_raw(:,1:2)-XY_center)*M_protein_RNA'+repmat(C_protein_RNA',size(XY_raw,1),1)*resolution_mismatch/resolution+XY_center;
        
    end
else
    xyz2 = xyz1;
end
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

