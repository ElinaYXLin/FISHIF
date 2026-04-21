function hist_fit_protein_diff(lf_name)
%clear all
close all

sub_folder = 'Histogram_protein_M/';
sub_folder2 = 'Histogram_protein_P/';
% sub_folder_nullo = 'Histogram_default/';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
data_tail = '_raw.xls';
data_tail2 = '_raw.xlsx';
output_add = '_spot_fit';
mask_area_add = '_mask_area';
mat_tail = '.mat';
fig_tail = '.fig';
hist_min = 0; 
hist_max = 4e5;
hist_bin = 5e2;
binw2 = 7.5e2;
fit_initial = [0.02,0.035,0.01,0.001,0.001,2.5e3,2e3,5e3,4e3,1];
fit_lower = [0,0,0,0,0,0,0,0,0,0];
fit_upper = [1,1,1,1,1,2.5e4,2.5e4,3.5e4,3e4,0.01];

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
        
    for list_J = 1:M1
        data_name = sub_list{list_J,3}(1:end-1);
        if exist([im_folder,sub_folder,data_name,data_tail],'file')
            raw_data = xlsread([im_folder,sub_folder,data_name,data_tail]);
        else
            raw_data = xlsread([im_folder,sub_folder,data_name,data_tail2]);
        end
        Inten_spot = prod(raw_data(:,1:3),2)*2*pi;
        n_hist = hist0(Inten_spot,x_hist_min,x_hist_max);
        
        load([im_folder,sub_folder,sub_list{list_J,3}(1:end-1),mask_area_add,mat_tail])
        mask_area1 = nnz(embryo_mask);
        
        if exist([im_folder,sub_folder2,data_name,data_tail],'file')
            raw_data2 = xlsread([im_folder,sub_folder2,data_name,data_tail]);
        else
            raw_data2 = xlsread([im_folder,sub_folder2,data_name,data_tail2]);
        end
        Inten_spot2 = prod(raw_data2(:,1:3),2)*2*pi;
        n_hist2 = hist0(Inten_spot2,x_hist_min,x_hist_max);
        
        load([im_folder,sub_folder2,sub_list{list_J,3}(1:end-1),mask_area_add,mat_tail])
        mask_area2 = nnz(embryo_mask);
        
%         n_hist = conv(n_hist,[0.5,1,0.5]/2,'same');
        x_hist = x_hist(1:end-1);
        y_hist = n_hist(1:end-1)/sum(n_hist(1:end-1));
        mgau = @(a0,a1,a2,a3,a4,b0,c0,b,c,d,x) a0*exp(-(x-b0).^2./2./c0.^2)+a1*exp(-(x-b).^2./2./c.^2)+a2*exp(-(x-2*b).^2./2./2./c.^2)+a3*exp(-(x-3*b).^2./3./2./c.^2)+a4*exp(-(x-4*b).^2./4./2./c.^2)+d;
        spot_fit = fit( x_hist',y_hist',mgau,'StartPoint',fit_initial,'Lower',fit_lower,'Upper',fit_upper);
        b = spot_fit.b;
        save([im_folder,sub_folder,data_name,output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
        figure
        bar(x_hist,y_hist,'b')
        hold on
        x_fit = [hist_min:(hist_max-hist_min)/1000:hist_max];
        y_fit = mgau(spot_fit.a0,spot_fit.a1,spot_fit.a2,spot_fit.a3,spot_fit.a4,spot_fit.b0,spot_fit.c0,spot_fit.b,spot_fit.c,spot_fit.d,x_fit);
        plot(x_fit,y_fit,'r')
        xlabel('Spot intensity (A.U.)')
        ylabel('Frequency')
        title(['Spot intensity histogram fit: ',data_name,', b = ',num2str(b)],'Interpreter','none')
        maximize(gcf)
        saveas(gcf,[im_folder,sub_folder,data_name,output_add,fig_tail])
        
%         b = 0;
%         if exist([im_folder,sub_folder_nullo]) ~= 7
%             mkdir([im_folder,sub_folder_nullo]);
%         end
%         save([im_folder,sub_folder_nullo,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
    end
end