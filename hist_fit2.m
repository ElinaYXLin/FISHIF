function hist_fit2(lf_name,out_fn)
%clear all
close all

sub_folder = 'Histogram_A/';
sub_folder_nullo = 'Histogram_default/';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
data_tail = '_raw.xls';
output_add = '_spot_fit';
mat_tail = '.mat';
fig_tail = '.fig';
hist_min = 0; 
hist_max = 4e5;
hist_bin = 2e3;
fit_initial = [0.035,0.001,0.0001,0.0001,1e4,1e4,1];
fit_lower = [0,0,0,0,0,0,0];

data_out = cell(0);
out_folder = 'FISH_calibration/';
out_tail = '.xls';
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
        data_name = [sub_list{list_J,3}(1:end-1),data_tail];
        raw_data = xlsread([im_folder,sub_folder,data_name]);
        Inten_spot = raw_data(:,1).*raw_data(:,2).*raw_data(:,3)*2*pi;
        [n_hist,x_hist] = hist(Inten_spot,[hist_min:hist_bin:(hist_max+hist_bin)]);
        x_hist = x_hist(1:end-1);
        y_hist = n_hist(1:end-1)/sum(n_hist(1:end-1));
        mgau = @(a1,a2,a3,a4,b,c,d,x) a1*exp(-(x-b).^2./2./c.^2)+a2*exp(-(x-2*b).^2./2./2./c.^2)+a3*exp(-(x-3*b).^2./3./2./c.^2)+a4*exp(-(x-4*b).^2./4./2./c.^2)+d;
        spot_fit = fit( x_hist',y_hist',mgau, 'StartPoint', fit_initial);
        b = spot_fit.b;
        c = spot_fit.c;
        
        data_out = cat(1,data_out,{folder_list{list_I,1},sub_list{list_J,3},b,c});
        
        save([im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
        figure
        bar(x_hist,y_hist,'b')
        hold on
        x_fit = [hist_min:(hist_max-hist_min)/1000:hist_max];
        y_fit = mgau(spot_fit.a1,spot_fit.a2,spot_fit.a3,spot_fit.a4,spot_fit.b,spot_fit.c,spot_fit.d,x_fit);
        plot(x_fit,y_fit,'r')
        xlabel('Spot intensity (A.U.)')
        ylabel('Frequency')
        title(['Spot intensity histogram fit: ',data_name,', b = ',num2str(b)],'Interpreter','none')
        saveas(gcf,[im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,fig_tail])
        
        b = 0;
        if exist([im_folder,sub_folder_nullo]) ~= 7
            mkdir([im_folder,sub_folder_nullo]);
        end
        save([im_folder,sub_folder_nullo,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(out_folder) ~= 7
    mkdir(out_folder);
end
xlswrite([out_folder,out_fn,out_tail],data_out)




