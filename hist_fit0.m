function hist_fit0(lf_name)
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
output_add = '_spot_fit';
mat_tail = '.mat';
fig_tail = '.fig';
hist_min = 0; 
hist_max = 4e5;
hist_bin = 2e3;
fit_initial = [0.035,0.001,0.0001,0.0001,1e4,1e4,1];
fit_lower = [0,0,0,0,0,0,0];

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
        b = 3600;
        if exist([im_folder,sub_folder]) ~= 7
            mkdir([im_folder,sub_folder]);
        end
        save([im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'b');
        
        b = 1;
        if exist([im_folder,sub_folder_nullo]) ~= 7
            mkdir([im_folder,sub_folder_nullo]);
        end
        save([im_folder,sub_folder_nullo,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'b');
    end
end