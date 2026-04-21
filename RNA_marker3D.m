function marker_RNA_profile = RNA_marker3D(RNA_fit,RNA2_fit,single_Inten,Inten_thresh,single_Inten2,Inten_thresh2,dxy_thresh,EL_info,image_folder,N_cycle,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to quantify # nascent mRNAs per loci using the fiducial marker information
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% marker_RNA_profile (Output): Loci expression data [EL_marker, #RNA, foci_ID, EL_RNA];
%% RNA_fit (input): RNA foci fitting result;
%% RNA2_fit (input): Fiducial marker foci fitting result;
%% single_Inten (input): Intensity of single mRNA (signal);
%% single_thresh (input): Threshold intensity of foci mRNA (signal);
%% single_Inten2 (input): Intensity of single mRNA (marker);
%% single_thresh2 (input): Threshold intensity of foci mRNA (marker);
%% dxy_thresh (input): Threshold of the distance between marker and foci;
%% EL_info (input): EL information of the embryo;
%% image_folder (input): Name of the image folder;
%% N_cycle (input): Cycle # of the embryo;
%% varargin (input): {flip_EL, resolutions, {EL_bin_min,EL_bin_max}}.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Input parameter extraction: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    flip_EL = false;
else
    flip_EL = varargin{1};
end

if isempty(varargin) || length(varargin) < 2
    rxyz = ones(1,3);
else
    rxyz = varargin{2};
end

if isempty(varargin) || length(varargin) < 3
    bin_min = 0:0.025:0.95;
    bin_max = 0.05:0.025:1;
else
    bin_min = varargin{3}{1};
    bin_max = varargin{3}{2};
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Foci data generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inten1 = 2*pi*prod(RNA_fit(:,1:3),2);
Inten2 = 2*pi*prod(RNA2_fit(:,1:3),2);

Itrue1 = Inten1 > Inten_thresh;
Itrue2 = Inten2 > Inten_thresh2;

NRNA1 = Inten1(Itrue1)/single_Inten;
NRNA2 = Inten2(Itrue2)/single_Inten2;

xyz1 = RNA_fit(Itrue1,6:8);
xyz2 = RNA2_fit(Itrue2,6:8);

%%% Calculate the EL of foci:
x0 = EL_info(1);
y0 = EL_info(2);
x1 = EL_info(3);
y1 = EL_info(4);
L2_extreme = EL_info(5);
EL1 = 1-dot((xyz1(:,[2,1])-repmat([x0,y0],size(xyz1,1),1)),repmat([(x1-x0),(y1-y0)],size(xyz1,1),1),2)/L2_extreme;
EL2 = 1-dot((xyz2(:,[2,1])-repmat([x0,y0],size(xyz2,1),1)),repmat([(x1-x0),(y1-y0)],size(xyz2,1),1),2)/L2_extreme;
if flip_EL
    EL1 = 1-EL1;
    EL2 = 1-EL2;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Foci matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find matching:
dxyz = pdist2(xyz1.*repmat(rxyz,size(xyz1,1),1),xyz2.*repmat(rxyz,size(xyz2,1),1));   %%% Distance between marker and foci
[dxyz0,I_foci] = min(dxyz);

%%% Assign # nacent mRNAs:
marker_RNA_profile = [EL2,NRNA1(I_foci'),I_foci',EL1(I_foci')];
marker_RNA_profile(dxyz0 > dxy_thresh,2) = 0;   %%% Exclude foci that are far away
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xbin,ybin,xerr,yerr,N_in] = equal_dist(marker_RNA_profile(:,1),marker_RNA_profile(:,2),bin_min,bin_max);
xbin = (bin_min+bin_max)/2;
figure(1)
    plot(marker_RNA_profile(:,1),marker_RNA_profile(:,2),'k.','DisplayName','Data')
    hold on
    errorbar(xbin,ybin,yerr,'Marker','none','Color','r','DisplayName','Bin')
    xlabel('EL')
    ylabel('# RNA / Marker')
    title(['TX quantification using fiducial marker: ',image_folder,', Cycle ',num2str(N_cycle)],'Interpreter','none')
    hold off






