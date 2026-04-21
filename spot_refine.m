function [vx,vy] = spot_refine(Inten,r,Ibin,rbin,image_folder,CH_name,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to display 3D distribution of spots for peak and r and %%%%%%
%% manually select a region for real spots. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vx (output): vertex x coordinates for ROI;
%% vy (output): vertex y coordinates for ROI;
%% Inten (input): peak intensity of spots;
%% r (input): radius of spots;
%% Ibin (input): intensity bin;
%% rbin (input): radius bin;
%% image_folder (input): image folder name;
%% CH_name (input): image channel name.
%% varargin (input): {Ibin width, rbin width};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2D distribution plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    n2D = hist3([Inten,r],{Ibin,rbin});
else
    n2D = hist3w([Inten,r],{Ibin,rbin},[varargin{1},varargin{2}]);
end
figure(1)
maximize(1)
    imagesc(Ibin,rbin,n2D')
    axis xy
    xlabel('Intensity')
    ylabel('Radius')
    title([image_folder,', Channel ',num2str(CH_name),': please select ROI'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ROI selection: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vx = zeros(0);
    vy = zeros(0);
    while ~~isempty(vx)
        [~,vx,vy] = roipoly;
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Display ROI: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    patch(vx,vy,ones(size(vx)),'FaceColor','none','EdgeColor','k','LineWidth',1)
