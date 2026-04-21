function beta_all = nucleus_protein_fit(lf_name,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to fit protein profile: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
fit_folder_protein2 = 'Histogram_protein_A1/';
fit_add = '_spot_fit';

bin_min = 0:0.025:0.95;   %%% bin_min for EL
bin_max = 0.05:0.025:1;   %%% bin_max for EL
xline0 = 0:0.01:1;
EL_range2 = [0,0.75];
expf = @(beta,x) beta(3)*exp(-(x-beta(2))./beta(1))./(1+exp(-(x-beta(2))./beta(1)))+beta(4);
% sigmaz = 1.3;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

beta_all = cell(0);

[N1,N2] = size(folder_list);
for list_I = 1:N1
    [sub_num, sub_list, sub_all] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
%     channel_name = eval(folder_list{list_I,5});
    if ~isempty(varargin) && ~isempty(varargin{1})
        list_J_all = varargin{1}{list_I};
    else
        list_J_all = 1:M1;
    end
    if size(sub_all,2) >= 17
        I_se = ~strcmp('F',sub_all(:,end));
%         I_se = strcmp('T',sub_all(:,end));
        list_J_all = list_J_all(I_se(list_J_all));
    end
    
    for list_J = list_J_all%1:M1
        load([folder_list{list_I,1},out_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_protein_profile_ab','resolution','resolutionz','quanti_p');
% %         
% %         load([folder_list{list_I,1},fit_folder_protein,sub_list{list_J,3}(1:end-1),fit_add,mat_tail],'b');   %%% load spot intensity fitting result
% %         Inten_protein0 = resolution^2*resolutionz*6.02e8*sqrt(2*pi)*sigmaz*b;
% % % %         b0 = b;
% % % %         load([folder_list{list_I,1},fit_folder_protein2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail],'b');   %%% load spot intensity fitting result
% % % %         b1 = b;
        
        EL0 = nucleus_protein_profile_ab(:,1);
        C0 = nucleus_protein_profile_ab(:,5)/1e-9;
% %         C0 = nucleus_protein_profile_ab(:,2);
% %         C1 = nucleus_protein_profile_ab(:,2)*quanti_p(1)/Inten_protein0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [xbin,ybin,xerr,yerr,N_in] = equal_dist(EL0,C0,bin_min,bin_max);
        Itrue = (xbin >= EL_range2(1)) & (xbin <= EL_range2(2)) & ~isnan(ybin);
        xbin0 = xbin(Itrue);
        ybin0 = ybin(Itrue);
% %         if isempty(protein_prefit) || size(protein_prefit,1) < i00 || protein_prefit(i00,1) <= 0 || all(protein_prefit(i00,2:end) == 0) || any(isnan(protein_prefit(i00,2:end)))
            beta0 = [-0.2,0.4,max(ybin0),0];
            try
                [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,expf,beta0);
            catch err
                beta1 = nan(size(beta0));
            end
% %         else
% %             beta1 = protein_prefit(i00,2:end);
% %         end
        beta_all = cat(1,beta_all,cat(2,{[folder_list{list_I,1},out_folder,sub_list{list_J,3}]},num2cell([beta1,max(ybin0)])));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot protein profile: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
            plot(EL0,C0,'.','DisplayName','Data')
% %             plot(EL0,C0,'.','DisplayName','Fluctuation method')
            hold on
% %             plot(EL0,C1/1e-9,'.','DisplayName','Spot method')
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Bin');
            plot(xline0,expf(beta1,xline0),'r-','DisplayName','Fit');
            text(0.1,40,['EL1 = ',num2str(beta1(1),'%5.2f'),char(10),'EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
            xlabel('EL')
            ylabel('[Protein] (nM)')
            title([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'Interpreter','none')
            legend('show')
    end
end


