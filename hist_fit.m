function hist_fit(lf_name,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to fit the RNA spot histogram and extract the single RNA intensity %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lf_name (input): Data folder for processing, can be in two formats: (a) input as a character array (eg. 'Duallist.xls'): the name of an Excel file containing the list of data folders;
%%                                                                     (b) input as a cell array (eg. {'111/';'222'}): manually list data folders (the delimiter must be ";", not ",");
%% varargin (input): 1. Name of the sub_folder for processing (default: 'Histogram_A/' (for RNA channel 1). For RNA channel 2, please use 'Histogram_A_RNA2/');
%%                   2. Type of processing (default (0): plot histogram and extract the single RNA intensity; 
%%                                                    1: replot threshold value for RNA spot recognition;
%%                                                    2: replot the histogram of RNA spot intensity)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
close all

% sub_folder = 'Histogram_A/';
% sub_folder = 'Histogram_A_gal4/';
sub_folder = 'Histogram_A_RNA2/';
r_th = 0;
sub_folder_nullo = 'Histogram_default/';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
data_tail = '_raw.xls';
show_tail = '_foci_th.fig';
output_add = '_spot_fit';
mat_tail = '.mat';
fig_tail = '.fig';
hist0_min = 0;
hist0_max = 5e3;
hist0_bin = 1e1;
fit0_initial = [1,100,50,0.05,500,300];
fit0_lower = [0,0,0,0,300,0];
fit0_higher = [1,500,150,0.5,2000,2000];

hist_min = 0; 
hist_max = 4e5;
hist_bin = 3e3;
% % fit_initial = [0.035,0.001,0.0001,0.0001,3e4,3e4,0,0.01,8e3,8e3];
% % fit_lower = [0,0,0,0,0,0,0,0,0,0];
% % fit_higher = [0.5,0.5,0.5,0.5,5e4,5e4,1,0.5,1.2e4,1.2e4];
fit_initial = [0.035,0.001,0.0001,0.0001,2.5e4,2.5e4,0,0.01,5e3,5e3];
fit_lower = [0,0,0,0,0,0,0,0,0,0];
fit_higher = [0.5,0.5,0.5,0.5,4e4,4e4,1,0.1,1e4,1e4];

% x_thresh0 = 175;

if ~isempty(varargin) && ~isempty(varargin{1})
    sub_folder = varargin{1};
end

if isempty(varargin) || length(varargin) < 2 || isempty(varargin{2})
    is_show = 0;    
else
    is_show = varargin{2};
end

%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name = lf_name;
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    im_folder = folder_list{list_I,1};
        
    for list_J = 1:M1
    if ~is_show
        data_name = [sub_list{list_J,3}(1:end-1),data_tail];
        raw_data = xlsread([im_folder,sub_folder,data_name]);
        raw_data = raw_data(sqrt(raw_data(:,2).*raw_data(:,3)) >= r_th,:);
        
        xxx = raw_data(:,1); yyy = sqrt(raw_data(:,2).*raw_data(:,3));
        [nbin,xbin] = hist(xxx,[hist0_min:hist0_bin:(hist0_max+hist0_bin)]);
        nbin = nbin(1:end-1);
        xbin = xbin(1:end-1);
        [~,Imax] = max(nbin);
        fit0_initial(2) = xbin(Imax); fit0_initial(5) = 2.5*xbin(Imax);
        fit0_higher(2) = 2*xbin(Imax); fit0_higher(5) = 5*xbin(Imax);
        fit0_lower(5) = xbin(Imax);
        f = fit(xbin',nbin'/max(nbin),'gauss2', 'StartPoint', fit0_initial, 'Upper', fit0_higher, 'Lower', fit0_lower);
        x0 = [f.b1,f.b2];
        c0 = [f.c1,f.c2]/sqrt(2);
        a0 = [f.a1,f.a2];
% % %         [~,Imax] = max(x0);
% % %         x_thresh = x0(Imax)-c0(Imax);
% %         [~,Imin] = min(x0);
% %         x_thresh = x0(Imin)+4*c0(Imin);
        dy = [abs(diff(feval(f,xbin))/(xbin(2)-xbin(1)));0];
        dyd = [abs(diff(feval(f,xbin))/(xbin(2)-xbin(1))-diff(feval(f,x0))/diff(x0));0];
        dy2 = [0;diff(diff(feval(f,xbin)));0];
        if max(dy(xbin > min(x0) & xbin <= max(x0) & dy2' > 0)) < 0
            [~,Imin2] = min(dy(xbin > min(x0) & xbin <= max(x0) & dy2' > 0));
        else
            [~,Imin2] = min(dyd(xbin > min(x0) & xbin <= max(x0) & dy2' > 0));
        end
        xbin2 = xbin(xbin > min(x0) & xbin <= max(x0) & dy2' > 0);
        
        if ~exist('x_thresh0') || isempty(x_thresh0)
            x_thresh = xbin2(Imin2);
        else
            x_thresh = x_thresh0;
        end
        
        y_thresh = mean(yyy)+2*std(yyy);
        
        Itrue = xxx > x_thresh;% & yyy < y_thresh;
        
% % %         [idx,c0] = kmeans(xxx,3);
% % %         [~,idx0] = min(c0(:,1));
% % % % %         [idy,c0] = kmeans(yyy,2);
% % % % %         [~,idy0] = sort(c0(:,1));
% % %         Itrue = idx ~= idx0;% & idy == idy0(1);
        Inten_spot = raw_data(Itrue,1).*raw_data(Itrue,2).*raw_data(Itrue,3)*2*pi;
        
        figure;
        subplot(2,1,1)
        plot(xxx,yyy,'k.','DisplayName','All spots'); 
        hold on; 
        plot(xxx(Itrue),yyy(Itrue),'r.','DisplayName','Selected spots');
% % %         hist3([xxx,yyy],'CdataMode','auto','EdgeColor','none','Ctrs',{0:2e1:1e5,0:0.02:10})
% % %         colorbar;xlim([0,4e3]);ylim([0,8]);
% % %         view(2)        
        xlabel('Spot intensity (A.U.)')
        ylabel('Spot radius (pixel)')
        legend('show')
        title('Spot gating','Interpreter','none')
        
        [n_hist,x_hist] = hist(Inten_spot,[hist_min:hist_bin:(hist_max+hist_bin)]);
        x_hist = x_hist(1:end-1);
        y_hist = n_hist(1:end-1)/sum(n_hist(1:end-1));
%         mgau = @(a1,a2,a3,a4,b,c,d,x) a1*exp(-(x-b).^2./2./c.^2)+a2*exp(-(x-2*b).^2./2./2./c.^2)+a3*exp(-(x-3*b).^2./3./2./c.^2)+a4*exp(-(x-4*b).^2./4./2./c.^2)+d;
        mgau = @(a1,a2,a3,a4,b,c,d,a0,b0,c0,x) a0*exp(-(x-b0).^2./2./c0.^2)+a1*exp(-(x-b).^2./2./c.^2)+a2*exp(-(x-2*b).^2./2./2./c.^2)+a3*exp(-(x-3*b).^2./3./2./c.^2)+a4*exp(-(x-4*b).^2./4./2./c.^2)+d;
% %         fit_initial = fit_initial0;
% %         [~,Imax] = max(y_hist);
% %         fit_initial(5:6) = x_hist(Imax); fit_initial(9:10) = x_hist(Imax)/2; 
% %         fit_higher = fit_higher0;
% %         fit_higher(5:6) = x_hist(Imax)*2; fit_higher(9:10) = x_hist(Imax); 
        
        spot_fit = fit( x_hist',y_hist',mgau, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower);
        b = spot_fit.b;
        save([im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
% %         figure;%maximize(gcf);
        subplot(2,1,2)
        bar(x_hist,y_hist,'b')
        hold on
        x_fit = [hist_min:(hist_max-hist_min)/1000:hist_max];
        y_fit = mgau(spot_fit.a1,spot_fit.a2,spot_fit.a3,spot_fit.a4,spot_fit.b,spot_fit.c,spot_fit.d,spot_fit.a0,spot_fit.b0,spot_fit.c0,x_fit);
        plot(x_fit,y_fit,'r')
        xlabel('Spot intensity (A.U.)')
        ylabel('Frequency')
        title(['Spot intensity histogram fit: ',data_name,', b = ',num2str(b)],'Interpreter','none')
        maximize(gcf)
        saveas(gcf,[im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,fig_tail])
        
        b = 0;
        if exist([im_folder,sub_folder_nullo]) ~= 7
            mkdir([im_folder,sub_folder_nullo]);
        end
        save([im_folder,sub_folder_nullo,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');

    elseif is_show == 1
        th_name = [sub_list{list_J,3}(1:end-1),show_tail];
        h0 = openfig([im_folder,sub_folder,th_name]);
%         fit_name = [sub_list{list_J,3}(1:end-1),output_add,fig_tail];
%         h0 = openfig([im_folder,sub_folder,fit_name]);
        maximize(h0)
        
        data_name = [sub_list{list_J,3}(1:end-1),data_tail];
        raw_data = xlsread([im_folder,sub_folder,data_name]);
        [sub_list{list_J,3}(1:end-1),', N = ',num2str(size(raw_data,1))]
        
    elseif is_show == 2
%         th_name = [sub_list{list_J,3}(1:end-1),show_tail];
%         h0 = openfig([im_folder,sub_folder,th_name]);
        fit_name = [sub_list{list_J,3}(1:end-1),output_add,fig_tail];
        h0 = openfig([im_folder,sub_folder,fit_name]);
        maximize(h0)
        
        data_name = [sub_list{list_J,3}(1:end-1),data_tail];
        raw_data = xlsread([im_folder,sub_folder,data_name]);
        [sub_list{list_J,3}(1:end-1),', N = ',num2str(size(raw_data,1))]
        
    end
    end
end

