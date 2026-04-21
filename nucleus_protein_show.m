function nucleus_protein_show(lf_name,varargin)
%clear all
close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
out_folder = 'Results/';
input_name = 'matchlist.xls';
mat_tail = '.mat';

fit_folder_protein = 'Histogram_protein_A/';
fit_add = '_spot_fit';

sigmaz = 1.3;
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
        load([folder_list{list_I,1},out_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_protein_profile_ab','resolution','resolutionz','quanti_p');
        
        load([folder_list{list_I,1},fit_folder_protein,sub_list{list_J,3}(1:end-1),fit_add,mat_tail],'b');   %%% load spot intensity fitting result
        Inten_protein0 = resolution^2*resolutionz*6.02e8*sqrt(2*pi)*sigmaz*b;
        
        EL0 = nucleus_protein_profile_ab(:,1);
        C0 = nucleus_protein_profile_ab(:,2);
        C1 = nucleus_protein_profile_ab(:,2)*quanti_p(1)/Inten_protein0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot protein profile: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
            plot(EL0,C0/1e-9,'.','DisplayName','Fluctuation method')
            hold on
            plot(EL0,C1/1e-9,'.','DisplayName','Spot method')
            xlabel('EL')
            ylabel('[Protein] (nM)')
            title([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
            legend('show')
    end
end


