function foci_RNA_show(lf_name,varargin)
%clear all
% close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
out_folder = 'Results/';
input_name = 'matchlist.xls';
mat_tail = '.mat';
hist_folder = 'Histogram_alignment/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_alignment_RNA2/';
fit_folder2 = 'Histogram_A_RNA2/';
fit_folder_protein = 'Histogram_protein_A/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';

bin_min = 0:0.025:0.95;
bin_max = 0.05:0.025:1;
ylim0 = [0,60];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
%     channel_name = eval(folder_list{list_I,5});
    if ~isempty(varargin) && ~isempty(varargin{1})
        list_J_all = varargin{1}{list_I};
    else
        list_J_all = 1:M1;
    end
    
    for list_J = list_J_all%1:M1
        Inten0 = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail],'b');   %%% load spot intensity fitting result
        I0 = prod(Inten0(:,1:3),2)*2*pi/b;
        
        load([folder_list{list_I,1},out_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'em_mask');
        EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask
        foci_xy = Inten0(:,[7,6]);
        x0 = EL_info(1);
        y0 = EL_info(2);
        x1 = EL_info(3);
        y1 = EL_info(4);
        L2_extreme = EL_info(5);

        EL0 = 1-dot((foci_xy-repmat([x0,y0],size(foci_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(foci_xy,1),1),2)/L2_extreme;
        
        [xbin,ybin,~,~] = equal_dist(EL0,I0,bin_min,bin_max);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot protein profile: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
            plot(EL0,I0,'.')
            hold on
            plot(xbin,ybin,'r')
            xlabel('EL')
            ylabel('Nascent RNA (#)')
            title(sub_list{list_J,3})
            ylim(ylim0)
    end
end


