function peak_parameter = ptrack_par(imstack0,imstack,mask_out)
tic
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% A function to track/fit local maxima to gaussian shapes: %%%%%%%
%% peak_parameter: fitting parameter values for each local maxima with the 
%%                 form: (I0,sigma_x,sigma_y,sigma_z,A,x0,y0,z0,theta,residual)
%% imstack0: Raw image stack
%% imstack: Prefiltered image stack
%% mask_out: 3D local maxima from the primary filtering step
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global all_offset xdata ydata xdata1 ydata1 raw_parameter isfit gau3 X0 Y0 Z0 xf fit_range Ipeak Nfit Jfit

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xrange = 6;  %%% Peak fitting range on x direction
yrange = 6;  %%% Peak fitting range on y direction
zrange = 0;  %%% Peak fitting range on z direction
sigmax0 = 2;   %%% Initial value of sigma_x
sigmay0 = 2;   %%% Initial value of sigma_y
sigmaz0 = 1;   %%% Initial value of sigma_z
crossxy = 0;   %%% Initial value of cross coefficient for "xy" terms
gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator
gau1_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,3)-x(',num2str(n),',3)).^2))+x(',num2str(n),',4)'];   %%% 1D single-gaussian model function text generator
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Coordinate matrix generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim0 = size(imstack0);
pdim0 = prod(dim0);
%[XC,YC,ZC] = ndgrid(1:dim0(1),1:dim0(2),1:dim0(3));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak sorting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_out(:,:,[1,end]) = false;
Ipeak = find(imstack.*mask_out);
[~,Itr] = sort(imstack(Ipeak),'descend');
Ipeak = Ipeak(Itr);
[X0,Y0,Z0] = ind2sub(dim0,Ipeak);   %%% Peak coordinates list
% X0 = XC(Ipeak);   %%% Peak coordinates x list
% Y0 = YC(Ipeak);   %%% Peak coordinates y list
% Z0 = ZC(Ipeak);   %%% Peak coordinates z list
vpeak = imstack0(Ipeak);   %%% Peak value list
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_parameter = zeros(length(Ipeak),10);
exstatus = false(length(Ipeak),1);
all_offset = zeros(pdim0,1);
% lb = [-1e9,-1e9,-1e9,-1e9,-1e9,-1e9,-1e9];
% ub = [1e9,1e9,1e9,1e9,1e9,1e9,1e9];

lb = [-inf,-inf,-inf,-inf,-inf,-inf,-inf];
ub = [inf,inf,inf,inf,inf,inf,inf];

options = optimset('Display','off');

%%% Arrange fitting groups: %%% ===========================================
% Dx = pdist(X0);   %%%%% distance matrix on x dimension
% Dy = pdist(Y0);   %%%%% distance matrix on y dimension
% Dz = pdist(Z0);   %%%%% distance matrix on z dimension
% Dxyz = squareform((pdist(X0) <= 2*xrange)&(pdist(Y0) <= 2*yrange)&(pdist(Z0) <= 2*zrange)); %%%%% neighbor matrix
isfit = false(size(X0));   %%% Peak fitting status
If_peak = cell(0);
fmax = 1;

for Ifit = 1:length(Ipeak)
    if ~isfit(Ifit)   %%% check whether the peak has been fitted
        Jfit = Ifit;   %%% instant peak in process
        Nfit = [];   %%% total peaks in fitting
        %%%%% Fitting range/peak arrangement: %%%%% -----------------------
        while Jfit
            Nfit = union(Nfit,Jfit);
            Dxyz = (pdist2(X0(Jfit),X0) <= 2*xrange)&(pdist2(Y0(Jfit),Y0) <= 2*yrange)&(pdist2(Z0(Jfit),Z0) <= 2*zrange);
            [~,Jfit] = find(Dxyz);
            Jfit = setdiff(Jfit,Nfit);
        end
        %%%%%% ------------------------------------------------------------
        isfit(Nfit) = true;
        If_peak = cat(1,If_peak,{Nfit});
        fmax = max(fmax,length(Nfit));
    end
end
clear Dxyz
%%% =======================================================================

        
%%% Fitting preparation: %%% ==============================================
%%%%% 2D multi-gaussian fit function generation:
gau2 = cell(1,fmax);
for I_gau = 1:fmax
    textfun = '@(x,xdata) ';
    for n = 1:I_gau
        textfun = [textfun,gau2_gen(n),'+'];
    end
    textfun = [textfun(1:(end-1)),';'];
    gau2{I_gau} = eval(textfun);
end
        
%%%%% Data points collection: 
raw_cell = cell(length(If_peak),1);
ex_cell = cell(length(If_peak),1);

t_status = true;
while t_status
    t_status = false;
    try
        matlabpool(6)
    catch
        t_status = true;
    end
end

parfor Ifit = 1:length(If_peak)
    Nfit = If_peak{Ifit};
    Nfit = reshape(Nfit,numel(Nfit),1);
    xdata = zeros(0,3);
    ydata = zeros(0,3);
    for I_temp = 1:length(Nfit)
        [Xtemp,Ytemp,Ztemp] = ndgrid(max((X0(Nfit(I_temp))-xrange),1):min((X0(Nfit(I_temp))+xrange),dim0(1)),max((Y0(Nfit(I_temp))-yrange),1):min((Y0(Nfit(I_temp))+yrange),dim0(2)),max((Z0(Nfit(I_temp))-zrange),1):min((Z0(Nfit(I_temp))+zrange),dim0(3)));
        xdata = union(xdata,[Xtemp(:),Ytemp(:),Ztemp(:)],'rows');
    end
    ydata =imstack0(sub2ind(dim0,xdata(:,1),xdata(:,2),xdata(:,3)));
%%%%% Fitting parameter initialization:
    xf0 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit)), X0(Nfit), 1./2./sigmay0.^2.*ones(size(Nfit)), Y0(Nfit), crossxy.*ones(size(Nfit)), min(ydata)*ones(size(Nfit))];
    xf10 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit)), Z0(Nfit), min(ydata)*ones(size(Nfit))];
%%% =======================================================================

%%% Multi-Gaussian fitting: %%% ===========================================
    [xf,resnorm,~,exitflag,~] = lsqcurvefit(gau2{length(Nfit)},xf0,xdata,ydata,lb,ub,options);
    xf1 = xf10;
    exitflag1 = 1;
    raw_cell{Ifit} = [xf(:,1:6),xf1(:,2:3),xf(:,7),sqrt(resnorm./length(ydata)*ones(size(Nfit)))];
    ex_cell{Ifit} = (exitflag > 0) & (exitflag1 > 0) & ((xf(:,2)+xf(:,4)) >= 0) & (4*xf(:,2).*xf(:,4) >= xf(:,6).*xf(:,6)) & (xf(:,3) > 0) & (xf(:,3) < dim0(1)) & (xf(:,5) > 0) & (xf(:,5) < dim0(2)) & (xf1(:,3) > 0) & (xf1(:,3) < dim0(3));   %%% exit status
%%% =======================================================================
end

matlabpool close

%%% Repackage data: %%% ===================================================
raw_parameter = cell2mat(raw_cell);
exstatus = cell2mat(ex_cell);
%%% =======================================================================
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fitting parameter output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_parameter = raw_parameter(exstatus,:);
peak_parameter(:,[1,5:8,10]) = raw_parameter(:,[1,9,3,5,8,10]);
peak_parameter(:,2) = 1./sqrt(raw_parameter(:,2)+raw_parameter(:,4)+sqrt((raw_parameter(:,2)-raw_parameter(:,4)).^2+raw_parameter(:,6).^2));
peak_parameter(:,3) = 1./sqrt(raw_parameter(:,2)+raw_parameter(:,4)-sqrt((raw_parameter(:,2)-raw_parameter(:,4)).^2+raw_parameter(:,6).^2));
peak_parameter(:,4) = 1./sqrt(2*raw_parameter(:,7));
peak_parameter(:,9) = atan2(raw_parameter(:,6),(raw_parameter(:,2)-raw_parameter(:,4)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

