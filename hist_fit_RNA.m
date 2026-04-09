function hist_fit_RNA(lf_name,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to fit the RNA spot histogram and extract the single RNA intensity %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lf_name (input): Data folder for processing, can be in two formats: (a) input as a character array (eg. 'Duallist.xls'): the name of an Excel file containing the list of data folders;
%%                                                                     (b) input as a cell array (eg. {'111/';'222'}): manually list data folders (the delimiter must be ";", not ",");
%% varargin (input): 1. Name of the sub_folder for processing (default: {'Histogram_A/','Histogram_P/'} (for RNA channel 1). For RNA channel 2, please use {'Histogram_A_RNA2/';'Histogram_P_RNA2/'});
%%                   2. Type of processing (default (0): plot histogram and extract the single RNA intensity; 
%%                                                    1: replot threshold value for RNA spot recognition;
%%                                                    2: replot the histogram of RNA spot intensity)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
close all

sub_folder = 'Histogram_A/';
sub_folder2 = 'Histogram_P/';
% sub_folder = 'Histogram_A_RNA2/';
% sub_folder2 = 'Histogram_P_RNA2/';
% sub_folder = 'Histogram_A_gal4/';
r_th = 0;
% % % sub_folder_nullo = 'Histogram_default/';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
data_tail = '_raw.xls';
data_tail2 = '_raw.xlsx';
show_tail = '_foci_th.fig';
output_add = '_spot_fit_new2';
mat_tail = '.mat';
fig_tail = '.fig';
hist0_min = 0;
hist0_max = 5e3;
hist0_bin = 1e1;
fit0_initial = [1,100,50,0.05,500,300];
fit0_lower = [0,0,0,0,300,0];
fit0_higher = [1,500,150,0.5,2000,2000];

hist_min = 0; 
hist_max = 2e5;
hist_bin = 1.5e3;
binw2 = 0.6e3;

% % fit_initial = [0.035,0.001,0.0001,0.0001,3e4,3e4,0,0.01,8e3,8e3];
% % fit_lower = [0,0,0,0,0,0,0,0,0,0];
% % fit_higher = [0.5,0.5,0.5,0.5,5e4,5e4,1,0.5,1.2e4,1.2e4];
% fit_initial = [0.035,0.01,0.005,0.001,2.5e4,2e4,0,0.01,0.7e4,0.5e4];
% fit_lower = [0,0,0,0,0,0,0,0,0,0];
% fit_higher = [0.5,0.5,0.5,0.5,4.5e4,4.5e4,1,0.1,1.5e4,1e4];
% fit_initial = [0.035,0.01,0.005,0.001,0.000,2e4,2e4,0,0.01,1e4,1e4];
% fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
% fit_higher = [1,1,1,1,0,0.4e5,0.4e5,1,1,0.6e4,0.6e4];
% fit_range = [0,2e4];
%% hb kr
fit_initial = [0.035,0.01,0.005,0.001,0.001,5e4,5e4,0,0.1,2e4,2e4];
fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
fit_higher = [0.5,0.5,0.5,0.5,0.5,10e4,10e4,1,1,5e4,5e4];
fit_range = [0,2.5e5];
% %% kni hcr
% fit_initial = [0.035,0.01,0.005,0.001,0.001,10e4,10e4,0,0.1,5e4,5e4];
% fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
% fit_higher = [0.5,0.5,0.5,0.5,0.5,20e4,20e4,1,1,10e4,10e4];
% fit_range = [0,2e5];

xlim0 = [0,3.5e3];
ylim0 = [0,6.5];
bin3D = {0:0.6e2:3.5e3,0:0.1:6.5};
% x_thresh0 = 175;

if ~isempty(varargin) && ~isempty(varargin{1})
    sub_folder = varargin{1}{1};
    sub_folder2 = varargin{1}{2};
    
    if strcmp(sub_folder,'Histogram_A_RNA2/')
%         hist_min = 0; 
%         hist_max = 5e4;
%         hist_bin = 0.3e3;
%         binw2 = 0.15e3;
% 
%         fit_initial = [0.035,0.001,0.001,0.001,0.001,0.5e4,0.5e4,0,0.01,0.3e4,0.3e4];
%         fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
%         fit_higher = [0.5,0.5,0.5,0.5,0.5,1.5e4,1.5e4,1,0.1,1e4,1e4];
%         fit_range = [0,5e4];
% 
%         xlim0 = [0,0.4e3];
%         ylim0 = [0,5];
%         bin3D = {0:1e1:0.4e3,0:0.1:10};
hist_min = 0; 
hist_max = 2.5e5;
hist_bin = 1.2e3;
binw2 = 0.6e3;

%% hb Kr
fit_initial = [0.035,0.01,0.005,0.001,0.001,8e4,8e4,0,0.1,4e4,4e4];
fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
fit_higher = [0.5,0.5,0.5,0.5,0.5,10e4,10e4,1,1,8e4,8e4];
fit_range = [0,2e5];



xlim0 = [0,1.5e4];
ylim0 = [0,6.5];
bin3D = {0:1e2:1.5e4,0:0.1:6.5};
    end
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
x_hist = hist_min:hist_bin:hist_max;
x_hist_min = x_hist-binw2;
x_hist_max = x_hist+binw2;

for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    im_folder = folder_list{list_I,1};
        
    for list_J = 2%1:M1
        if ~is_show
            figure
            maximize(gcf)
            a(1) = subplot(2,2,1);
            a(2) = subplot(2,2,2);
            a3 = subplot(2,2,3:4);

            t_back = false;
            
            data_name = sub_list{list_J,3}(1:end-1);
            if exist([im_folder,sub_folder,data_name,data_tail],'file')
                raw_data = xlsread([im_folder,sub_folder,data_name,data_tail]);
            else
                raw_data = xlsread([im_folder,sub_folder,data_name,data_tail2]);
            end

            if exist([im_folder,sub_folder2,data_name,data_tail],'file')
                raw_data2 = xlsread([im_folder,sub_folder2,data_name,data_tail]);
            elseif exist([im_folder,sub_folder2,data_name,data_tail2],'file')
                raw_data2 = xlsread([im_folder,sub_folder2,data_name,data_tail2]);
            else
                raw_data2 = zeros(1,10);
            end

% %             Itrue1 = raw_data(:,2)./raw_data(:,3) > 0.4;
% %             Itrue2 = raw_data2(:,2)./raw_data2(:,3) > 0.4;
% %             X = [raw_data(Itrue1,1),sqrt(raw_data(Itrue1,2).*raw_data(Itrue1,3))];
% %             X2 = [raw_data2(Itrue2,1),sqrt(raw_data2(Itrue2,2).*raw_data2(Itrue2,3))];
            
            X = [raw_data(:,1),sqrt(raw_data(:,2).*raw_data(:,3))];
            X2 = [raw_data2(:,1),sqrt(raw_data2(:,2).*raw_data2(:,3))];
            
    % % % %         X = [raw_data(:,1),raw_data(:,4)];
    % % % %         X2 = [raw_data2(:,1),raw_data2(:,4)];
    % %         X = [raw_data(:,1),raw_data(:,2)./raw_data(:,3)];
    % %         X2 = [raw_data2(:,1),raw_data2(:,2)./raw_data2(:,3)];

            axes(a(1))
    %             hist3(X,'CdataMode','auto','EdgeColor','none','Ctrs',{0:5e1:1e5,0:0.1:10})
                [N0,C0] = hist3(X,'CdataMode','auto','EdgeColor','none','Ctrs',bin3D);
                imagesc(C0{1},C0{2},log10(N0)');
                axis xy
                colorbar
    %             view(2)
% %                     hold on
    % %             plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2)
                xlim(xlim0)
                ylim(ylim0)
                xlabel('Peak intensity (A.U.)')
                ylabel('/sigma (pixel)')
                title('Anterior spots');
            axes(a(2))
    %             hist3(X2,'CdataMode','auto','EdgeColor','none','Ctrs',{0:5e1:1e5,0:0.1:10})
                [N0,C0] = hist3(X2,'CdataMode','auto','EdgeColor','none','Ctrs',bin3D);
                imagesc(C0{1},C0{2},log10(N0)');
                axis xy
                colorbar
    %             view(2)
% %                     hold on
    % %             plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2)
                xlim(xlim0)
                ylim(ylim0)
                xlabel('Peak intensity (A.U.)')
                ylabel('/sigma (pixel)')
                title('Posterior spots');
            linkaxes(a)

            h1 = [];
            while ~t_back
                cla(a3)
                t_gate0 = false;
                while ~t_gate0
                    if ~exist('Inten_lim') || isempty(Inten_lim)
                        axes(a(1))
                        xlim(xlim0)
                        ylim(ylim0)
                        xlabel('Peak intensity (A.U.)')
                        ylabel('/sigma (pixel)')
                        title('Please pick up the boundary');
        % %                 [x0,y0] = ginput(2);
        % %                 p0 = polyfit(y0,x0,1); 
%                         [x0,~] = ginput(1);
%                         p0(1) = 0; p0(2) = x0;
%                     else
%                         p0(1) = 0; p0(2) = Inten_lim;
%                     end
                    
                       [x0,y0] = ginput(2);
                       if length(x0) > 1
                          p0 = polyfit(y0,x0,1); 
                       else
                          p0(1) = 0; p0(2) = x0;
                       end
                    else
                     p0(1) = 0; p0(2) = Inten_lim;
                     t_gate0 = true;
                    end

                    delete(h1);
                    axes(a(1))
                        hold on
                        h1(1) = plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2);
            %             xlim(xlim0)
                    axes(a(2))
                        hold on
                        h1(2) = plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2);
                        xlim(xlim0)

                    if ~t_gate0
                        answer0 = questdlg('Is the gating OK','Gating confirmation','Yes','No','Yes');
                        t_gate0 = strcmpi('Yes',answer0);
                    end
                end


        %% Fit histogram: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
                Itrue = polyval(p0,X(:,2)) <= X(:,1);
                Inten_spot = prod(raw_data(Itrue,1:3),2)*2*pi;
                n_hist = hist0(Inten_spot,x_hist_min,x_hist_max);

        % %         Itrue2 = polyval(p0,X2(:,2)) <= X2(:,1);
        % %         Inten_spot2 = prod(raw_data2(Itrue2,1:3),2)*2*pi;
        % %         n_hist2 = hist0(Inten_spot2,x_hist_min,x_hist_max);

        % % %         Itrue2 = raw_data2(:,1) >= Inten_lim;
        % % %         Inten_spot2 = prod(raw_data2(Itrue2,1:3),2)*2*pi;
        % % %         n_hist2 = hist0(Inten_spot2,x_hist_min,x_hist_max);
        % % %         
        % % %         Itrue3 = raw_data3(:,1) >= Inten_lim;
        % % %         Inten_spot3 = prod(raw_data3(Itrue3,1:3),2)*2*pi;
        % % %         n_hist3 = hist0(Inten_spot3,x_hist_min,x_hist_max);

        %         n_hist = conv(n_hist,[0.5,1,0.5]/2,'same');
        % %         x_hist = x_hist(1:end-1);
        % %         y_hist = n_hist(1:end-1)/sum(n_hist(1:end-1));
        % % % %         y_hist = n_hist-n_hist2; y_hist = y_hist/sum(y_hist);
                y_hist = n_hist/sum(n_hist);%-n_hist2/sum(n_hist2)/4;
                mgau = @(a1,a2,a3,a4,a5,b,c,d,a0,b0,c0,x) a0*exp(-(x-b0).^2./2./c0.^2)+a1*exp(-(x-b).^2./2./c.^2)+a2*exp(-(x-2*b).^2./2./2./c.^2)+a3*exp(-(x-3*b).^2./3./2./c.^2)+a4*exp(-(x-4*b).^2./4./2./c.^2)+a5*exp(-(x-5*b).^2./5./2./c.^2)+d;
                Itrue0 = x_hist >= fit_range(1) & x_hist <= fit_range(2);
                spot_fit = fit( x_hist(Itrue0)',y_hist(Itrue0)',mgau, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower);
                b = spot_fit.b;

                axes(a3)
                    bar(x_hist,y_hist,'b')
                    hold on
                    x_fit = [hist_min:(hist_max-hist_min)/1000:hist_max];
                    y_fit = feval(spot_fit,x_fit);
                    plot(x_fit,y_fit,'r')
                    hold off
                    xlabel('Spot intensity (A.U.)')
                    ylabel('Frequency')
                    title(['Spot intensity histogram fit: ',data_name,', b = ',num2str(b)],'Interpreter','none')
            
                if ~t_back
                    answer0 = questdlg('Is the fitting OK','Fitting confirmation','Yes','No','Yes');
                    t_back = strcmpi('Yes',answer0);
                end

            end
            
            saveas(gcf,[im_folder,sub_folder,data_name,output_add,fig_tail])
            save([im_folder,sub_folder,data_name,output_add,mat_tail],'spot_fit','x_hist','y_hist','b');

    % % %         b = 0;
    % % %         if exist([im_folder,sub_folder_nullo]) ~= 7
    % % %             mkdir([im_folder,sub_folder_nullo]);
    % % %         end
    % % %         save([im_folder,sub_folder_nullo,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');

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

