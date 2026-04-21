function hist_gate_en(lf_name)
%clear all
close all

sub_folder = 'Histogram_en/';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
data_tail = '_raw.xls';
data_tail2 = '_raw.xlsx';
output_add = '_spot_fit';
mat_tail = '.mat';
fig_tail = '.fig';
hist_min = 0; 
hist_max = 4e5;
hist_bin = 2e3;
binw2 = 2e3;
fit_initial = [0.02,0.1,0.04,0.01,0.01,2e4,2e4,5e4,5e4,0];
fit_lower = [0,0,0,0,0,0,0,0,0,0];
fit_upper = [0.1,1,1,1,1,4e4,4e4,10e4,10e4,0];

xlim0 = [0,5e3];
ylim0 = [0,5];
% Inten_lim = 500;

%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name = lf_name;
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('cgal4new',folder_list(:,6)),:);
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
        
    for list_J = 1:M1
        data_name = sub_list{list_J,3}(1:end-1);
        if exist([im_folder,sub_folder,data_name,data_tail],'file')
            raw_data = xlsread([im_folder,sub_folder,data_name,data_tail]);
        else
            raw_data = xlsread([im_folder,sub_folder,data_name,data_tail2]);
        end
        
        X = [raw_data(:,1),sqrt(raw_data(:,2).*raw_data(:,3))];
        
        figure('Name',[im_folder,sub_folder,data_name])
        maximize(gcf)
        a(1) = subplot(1,2,1);
%             hist3(X,'CdataMode','auto','EdgeColor','none','Ctrs',{0:5e1:1e5,0:0.1:10})
            [N0,C0] = hist3(X,'CdataMode','auto','EdgeColor','none','Ctrs',{0:5e1:1e5,0:0.1:10});
            imagesc(C0{1},C0{2},log10(N0)');
            axis xy
            colorbar
%             view(2)
            hold on
%             plot3(Inten_lim*ones(size(ylim0)),ylim0,1000*ones(size(ylim0)),'w-','LineWidth',2)
            xlim(xlim0)
            ylim(ylim0)
            xlabel('Peak intensity (A.U.)')
            ylabel('/sigma (pixel)')
            title('Please pick up the boundary');
            
        t_gate0 = false;
        h1 = [];
        while ~t_gate0
            [x0,y0] = ginput(2);
            if length(x0) > 1
                p0 = polyfit(y0,x0,1); 
            else
                p0(1) = 0; p0(2) = x0;
            end
        
% %         axes(a(1))
            hold on
            delete(h1);
            h1 = plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2);
            xlim(xlim0)
            ylim(ylim0)
            
            answer0 = questdlg('Is the gating OK','Gating confirmation','Yes','No','Yes');
            t_gate0 = strcmpi('Yes',answer0);
        end
%% Fit histogram: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Itrue = polyval(p0,X(:,2)) <= X(:,1);
        Inten_spot = prod(raw_data(Itrue,1:3),2)*2*pi;
        n_hist = hist0(Inten_spot,x_hist_min,x_hist_max);
        
        y_hist = n_hist/sum(n_hist);%-n_hist2/sum(n_hist2)/4;
        
        save([im_folder,sub_folder,data_name,output_add,mat_tail],'p0');
        
        subplot(1,2,2)
            bar(x_hist,y_hist,'b')
            xlabel('Spot intensity (A.U.)')
            ylabel('Frequency')
            title('Spot intensity histogram','Interpreter','none')
        saveas(gcf,[im_folder,sub_folder,data_name,output_add,fig_tail])
    end
end