function [peak_parameter,varargout] = ptrack(imstack0,imstack,mask_out,varargin)
tic
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% A function to track/fit local maxima to gaussian shapes: %%%%%%%
%% %%%%%%% version: 11/20/2018: improve the speed of peak search         %% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% peak_parameter (Output): fitting parameter values for each local maxima (I0,sigma_x,sigma_y,sigma_z,A,x0,y0,z0,theta,residual);
%% varargout (Output): xf_all: raw fitting result from lsqcurvefit;
%% imstack0 (Input): Raw image stack;
%% imstack (Input): Prefiltered image stack;
%% mask_out (Input): 3D local maxima from the primary filtering step;
%% varargin (Input): (1) mask_out2D: 2D-only maxima;
%%                   (2) z00: boundary layers that do not need spot detection (default = 1);
%%                   (3) [xrange,yrange]: Peak fitting range on x and y directions (default = 5);
%%                   (4) Nfitmax: Maximum number of spot candidate for simultaneous fitting (default = 10);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global all_offset xdata ydata xdata1 ydata1 raw_parameter isfit gau3 X0 Y0 Z0 xf fit_range Ipeak Nfit Jfit

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) && ~isempty(varargin{1})
    mask_out2D = varargin{1};
else
    mask_out2D = false(size(mask_out));
end
if length(varargin) > 1 && ~isempty(varargin{2})   %%% boundary layers that do not need spot detection
    z00 = varargin{2};
else
    z00 = 1;
end
if length(varargin) > 2 && length(varargin{3}) > 1
    xrange = varargin{3}(1);  %%% Peak fitting range on x direction
    yrange = varargin{3}(2);  %%% Peak fitting range on y direction
    if length(varargin{3}) > 2
        zrange = varargin{3}(3);  %%% Peak fitting range on z direction
    else
        zrange = 2;
    end
elseif length(varargin) > 2 && length(varargin{3}) == 1
    xrange = varargin{3}(1);
    yrange = varargin{3}(1);
    zrange = 2;  %%% Peak fitting range on z direction
else
    xrange = 3;
    yrange = 3;
    zrange = 2;  %%% Peak fitting range on z direction
end
if length(varargin) >= 4 && ~isempty(varargin{4})
    Nfitmax = varargin{4};
else
    Nfitmax = 10;
end
sigmax0 = 2;   %%% Initial value of sigma_x
sigmay0 = 2;   %%% Initial value of sigma_y
sigmaz0 = 1;   %%% Initial value of sigma_z
crossxy = 0;   %%% Initial value of cross coefficient for "xy" terms
gau3_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2+x(',num2str(n),',8)*(xdata(:,3)-x(',num2str(n),',9)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator
gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator
gau1_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,3)-x(',num2str(n),',3)).^2))+x(',num2str(n),',4)'];   %%% 1D single-gaussian model function text generator
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Coordinate matrix generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear varargin
dim0 = size(imstack0);
pdim0 = prod(dim0);
%[XC,YC,ZC] = ndgrid(1:dim0(1),1:dim0(2),1:dim0(3));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak sorting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_out(:,:,[1:z00,end+1-z00:end]) = false;
Ipeak = find(mask_out);
[~,Itr] = sort(imstack(Ipeak),'descend');
Ipeak = Ipeak(Itr);
[X0,Y0,Z0] = ind2sub(dim0,Ipeak);   %%% Peak coordinates list
% X0 = XC(Ipeak);   %%% Peak coordinates x list
% Y0 = YC(Ipeak);   %%% Peak coordinates y list
% Z0 = ZC(Ipeak);   %%% Peak coordinates z list
vpeak = double(imstack0(Ipeak));   %%% Peak value list
isfit = false(size(X0));   %%% Peak fitting status

mask_out2D(:,:,[1,end]) = false;
Ipeak2D = find(mask_out2D);
[X02D,Y02D,Z02D] = ind2sub(dim0,Ipeak2D);   %%% Peak coordinates list
vpeak2D = double(imstack0(Ipeak2D));   %%% Peak value list
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_parameter = zeros(length(Ipeak),10);
xf_all = zeros(0);
exstatus = false(length(Ipeak),1);
% all_offset = zeros(pdim0,1);
all_offset = 0;
% lb = [-1e9,-1e9,-1e9,-1e9,-1e9,-1e9,-1e9];
% ub = [1e9,1e9,1e9,1e9,1e9,1e9,1e9];

% lb = [-inf,-inf,-inf,-inf,-inf,-inf,-inf];
% ub = [inf,inf,inf,inf,inf,inf,inf];
lb = [0,0,0,0,0,-20,0,0,0];
ub = [5e5,10,size(imstack,1),10,size(imstack,2),20,1e5,10,size(imstack,3)];

options = optimset('Display','off');
IIfit = 0;
num_peak = length(Ipeak);
for Ifit = 1:num_peak
%     Ifit
    if ~isfit(Ifit)   %%% check whether the peak has been fitted
        Jfit = Ifit;   %%% instant peak in process
        Nfit = [];   %%% total peaks in fitting
        Jfit2D = [];   %%% instant 2D peak in process
        Nfit2D = [];   %%% total 2D peaks in fitting
        %%% Fitting range/peak arrangement: %%%============================
%         fit_range = false(size(XC));
%         fit_range1 = false(size(XC));
% %         fit_range0 = zeros(0,3);
        while ~isempty(Jfit) && size([Nfit;Nfit2D],1) <= Nfitmax
            Nfit = cat(1,Nfit,Jfit);
            if (~isempty(Nfit2D)) && (~isempty(Jfit2D))
                Nfit2D = cat(1,Nfit2D,Jfit2D);
            end
            
            rx = pdist2(X0,X0(Nfit));
            ry = pdist2(Y0,Y0(Nfit));
            rz = pdist2(Z0,Z0(Nfit));
            Jpeak = find(any(rx <= 2*xrange & ry <= 2*yrange & rz <= 2*zrange,2));
            Jfit = setdiff(setdiff(Jpeak,Nfit),find(isfit));
            
            if ~isempty(X02D)
                rx2D = pdist2(X02D,X0(Nfit));
                ry2D = pdist2(Y02D,Y0(Nfit));
                rz2D = pdist2(Z02D,Z0(Nfit));
                Jpeak2D = find(any(rx2D <= 2*xrange & ry2D <= 2*yrange & rz2D <= 2*zrange,2));
                Jfit2D = setdiff(Jpeak2D,Nfit2D);
            else
                Jpeak2D = zeros(0);
                Jfit2D = zeros(0);
            end
        end
        Nfit2D = cat(1,Nfit2D,Jfit2D);
       
        %%% ===============================================================
        
        %%% Fitting preparation: %%%=======================================
        %%%%% 3D multi-gaussian fit function generation:
%         textfun = 'gau3 = @(x,xdata) ';
%         for n = 1:length(Nfit)
%             textfun = [textfun,gau3_gen(n),'+'];
%         end
        textfun3 = 'gau3 = @(x,xdata) ';
% %         textfun = 'gau2 = @(x,xdata) ';
% %         textfun1 = 'gau1 = @(x,xdata) ';
        for n = 1:(length(Nfit)+length(Nfit2D))
            textfun3 = [textfun3,gau3_gen(n),'+'];
% %             textfun = [textfun,gau2_gen(n),'+'];
% %             textfun1 = [textfun1,gau1_gen(n),'+'];
        end
        textfun3 = [textfun3(1:(end-1)),';'];
% %         textfun = [textfun(1:(end-1)),';'];
% %         textfun1 = [textfun1(1:(end-1)),';'];
        eval(textfun3);
% %         eval(textfun);
% %         eval(textfun1);
        %%%%% Data points collection: 
%         range1D = find(fit_range);
%         range11D = find(fit_range1);
%         xdata = [XC(range1D),YC(range1D),ZC(range1D)];   %%%%% Collect the coordinates of data points
%         ydata = imstack0(range1D)-all_offset(range1D);   %%%%% Collect the intensity values of data points
        xdata = zeros(0,3);
        ydata = zeros(0,3);
        for I_temp = 1:length(Nfit)
            [Xtemp,Ytemp,Ztemp] = ndgrid(max((X0(Nfit(I_temp))-xrange),1):min((X0(Nfit(I_temp))+xrange),dim0(1)),max((Y0(Nfit(I_temp))-yrange),1):min((Y0(Nfit(I_temp))+yrange),dim0(2)),max((Z0(Nfit(I_temp))-zrange),1):min((Z0(Nfit(I_temp))+zrange),dim0(3)));
            xdata = union(xdata,[Xtemp(:),Ytemp(:),Ztemp(:)],'rows');
        end
        for I_temp = 1:length(Nfit2D)
            [Xtemp,Ytemp,Ztemp] = ndgrid(max((X02D(Nfit2D(I_temp))-xrange),1):min((X02D(Nfit2D(I_temp))+xrange),dim0(1)),max((Y02D(Nfit2D(I_temp))-yrange),1):min((Y02D(Nfit2D(I_temp))+yrange),dim0(2)),max((Z02D(Nfit2D(I_temp))-zrange),1):min((Z02D(Nfit2D(I_temp))+zrange),dim0(3)));
            xdata = union(xdata,[Xtemp(:),Ytemp(:),Ztemp(:)],'rows');
        end
        ydata = double(imstack0(sub2ind(dim0,xdata(:,1),xdata(:,2),xdata(:,3))));
%         xdata1 = [XC(range11D),YC(range11D),ZC(range11D)];   %%%%% Collect the coordinates of data points
%         ydata1 = imstack(range11D)-all_offset(range11D);   %%%%% Collect the intensity values of data points

        %%%%% Fitting parameter initialization:
        xf0 = [vpeak(Nfit)-all_offset-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit)), X0(Nfit), 1./2./sigmay0.^2.*ones(size(Nfit)), Y0(Nfit), crossxy.*ones(size(Nfit)), min(ydata)*ones(size(Nfit)), 1./2./sigmaz0.^2.*ones(size(Nfit)), Z0(Nfit)];
% %         xf0 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit)), X0(Nfit), 1./2./sigmay0.^2.*ones(size(Nfit)), Y0(Nfit), crossxy.*ones(size(Nfit)), min(ydata)*ones(size(Nfit))];
% %         xf10 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit)), Z0(Nfit), min(ydata)*ones(size(Nfit))];
% % %         xf0 = [vpeak(Nfit)-all_offset-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit)), X0(Nfit), 1./2./sigmay0.^2.*ones(size(Nfit)), Y0(Nfit), crossxy.*ones(size(Nfit)), min(ydata)*ones(size(Nfit))];
% % %         xf10 = [vpeak(Nfit)-all_offset-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit)), Z0(Nfit), min(ydata)*ones(size(Nfit))];

        xf02D = [vpeak2D(Nfit2D)-all_offset-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit2D)), X02D(Nfit2D), 1./2./sigmay0.^2.*ones(size(Nfit2D)), Y02D(Nfit2D), crossxy.*ones(size(Nfit2D)), min(ydata)*ones(size(Nfit2D)), 1./2./sigmaz0.^2.*ones(size(Nfit2D)), Z02D(Nfit2D)];
% %         xf02D = [vpeak2D(Nfit2D)-all_offset(Ipeak2D(Nfit2D))-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit2D)), X02D(Nfit2D), 1./2./sigmay0.^2.*ones(size(Nfit2D)), Y02D(Nfit2D), crossxy.*ones(size(Nfit2D)), min(ydata)*ones(size(Nfit2D))];
% %         xf102D = [vpeak2D(Nfit2D)-all_offset(Ipeak2D(Nfit2D))-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit2D)), Z02D(Nfit2D), min(ydata)*ones(size(Nfit2D))];
% % %         xf02D = [vpeak2D(Nfit2D)-all_offset-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit2D)), X02D(Nfit2D), 1./2./sigmay0.^2.*ones(size(Nfit2D)), Y02D(Nfit2D), crossxy.*ones(size(Nfit2D)), min(ydata)*ones(size(Nfit2D))];
% % %         xf102D = [vpeak2D(Nfit2D)-all_offset-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit2D)), Z02D(Nfit2D), min(ydata)*ones(size(Nfit2D))];
        
        lb0 = repmat(lb,size(xf0,1),1);
        ub0 = repmat(ub,size(xf0,1),1);
        if ~isempty(Nfit)
            lb0(:,[3,5,9]) = [X0(Nfit)-xrange,Y0(Nfit)-yrange,Z0(Nfit)-zrange];
            ub0(:,[3,5,9]) = [X0(Nfit)+xrange,Y0(Nfit)+yrange,Z0(Nfit)+zrange];
        end
        
        lb02D = repmat(lb,size(xf02D,1),1);
        ub02D = repmat(ub,size(xf02D,1),1);
        if ~isempty(Nfit2D)
            lb02D(:,[3,5,9]) = [X02D(Nfit2D)-xrange,Y02D(Nfit2D)-yrange,Z02D(Nfit2D)-zrange];
            ub02D(:,[3,5,9]) = [X02D(Nfit2D)+xrange,Y02D(Nfit2D)+yrange,Z02D(Nfit2D)+zrange];
        end
        %%% ===============================================================
        
        if size([xf0;xf02D],1) <= Nfitmax && ~isempty([xf0;xf02D]) && ~isempty(xdata) && ~isempty(ydata)
            try
            %%% Multi-Gaussian fitting: %%%====================================
            [xf,resnorm,~,exitflag,~] = lsqcurvefit(gau3,[xf0;xf02D],xdata,ydata,[lb0;lb02D],[ub0;ub02D],options);
            %[xf,resnorm,~,exitflag,~] = lsqcurvefit(gau2,xf0,xdata,ydata);
            %[xf1,~,~,exitflag1,~] = lsqcurvefit(gau1,xf10,xdata1,ydata1);
            xf1 = xf(:,[1,8,9,7]);
    % %         xf1 = xf10;
            xf_all = cat(1,xf_all,xf);
            xf = xf(1:length(Nfit),:);
            xf(:,7) = sum(xf(:,7));
            xf1 = xf1(1:length(Nfit),:);
            exitflag1 = 1;
            raw_parameter(Nfit,:) = [xf(:,1:6),xf1(:,2:3),xf(:,7),sqrt(resnorm./length(ydata)*ones(size(Nfit)))];
            %raw_parameter(Nfit,:) = [xf(:,1:6),repmat(Z0(Ifit),size(xf,1),1),xf1(:,3),xf(:,7),sqrt(resnorm./length(ydata)*ones(size(Nfit)))];
            exstatus(Nfit) = (exitflag > 0) & (exitflag1 > 0) & ((xf(:,2)+xf(:,4)) >= 0) & (4*xf(:,2).*xf(:,4) >= xf(:,6).*xf(:,6)) & (xf(:,3) > 0) & (xf(:,3) < dim0(1)) & (xf(:,5) > 0) & (xf(:,5) < dim0(2)) & (xf(:,9) > 0) & (xf(:,9) < dim0(3));   %%% exit status
            %xf_offset = [xf(:,1:8),zeros(size(xf,1),1)];
            %all_offset = all_offset+gau3(xf_offset,[XC([1:pdim0]'),YC([1:pdim0]'),ZC([1:pdim0]')]);
            %%% ===============================================================
            catch
                [xdata,ydata]
                [xf0;xf02D]
            end
        else
            [xf0;xf02D]
        end
        isfit(Nfit) = true;
    end
    IIfit = IIfit+1;
    if IIfit >= 100
        disp(['Finishing: ',num2str(Ifit),'/',num2str(num_peak),', time: ',num2str(toc)])
        IIfit = 0;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fitting parameter output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_parameter = raw_parameter(exstatus,:);
peak_parameter(:,[1,5:8,10]) = raw_parameter(:,[1,9,3,5,8,10]);
peak_parameter(:,2) = 1./sqrt(raw_parameter(:,2)+raw_parameter(:,4)+sqrt((raw_parameter(:,2)-raw_parameter(:,4)).^2+raw_parameter(:,6).^2));
peak_parameter(:,3) = 1./sqrt(raw_parameter(:,2)+raw_parameter(:,4)-sqrt((raw_parameter(:,2)-raw_parameter(:,4)).^2+raw_parameter(:,6).^2));
peak_parameter(:,4) = 1./sqrt(2*raw_parameter(:,7));
peak_parameter(:,9) = atan2(raw_parameter(:,6),(raw_parameter(:,2)-raw_parameter(:,4)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout = {xf_all};
%clear XC YC ZC
toc
