clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to save the crosstalk parameters from analysis results %%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNA2list.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'Results/';
mat_tail = '.mat';

all_color0 = {'DAPI','A488','TMR','A594','A647'};

old_new = {'','','_new'};
I_old_new0 = 3;
corr_add = '_corr';
cross_add = '_cross';
mismatch_add = '_0mis';
intercept_add = '_0int';

save_folder = 'Calibration_crosstalk/';
% % save_name = '02072016';
save_name = '01052020';

list_J0 = 2:3;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extract crosstalk analysis results: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P0 = cell(0);
Cin = cell(0);
Cout = cell(0);

for list_I = 1:N1
    all_foci = cell(0);
    all_foci_xy = cell(0);
    all_single = cell(0);
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    channel_name = eval(folder_list{list_I,5});
    [M1,M2] = size(sub_list);
    
    for list_J = list_J0%1:M1
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,I_old_new0}(1:(length(sub_list{list_J,I_old_new0})-1)),mismatch_add,intercept_add,mat_tail]);
        
        for ii = 1:length(channel_in)
            for jj = 1:length(channel_out)
                Cin = cat(1,Cin,all_color(channel_in(ii)));
                Cout = cat(1,Cout,all_color(channel_out(jj)));
                P0 = cat(1,P0,p_all(ii,jj));
            end
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Construct the crosstalk matrices and output: %%%%%%%%%%%%%%%%%%%%%%%%%%%
CrossA = diag(ones(size(all_color0)));
CrossB = zeros(size(all_color0))';

for ii = 1:length(all_color0)
    Iin = strcmp(all_color0{ii},Cin);
    
    for jj = 1:length(all_color0)
        Iout = strcmp(all_color0{jj},Cout);
        
        if any(Iin & Iout)
            p0_mean = mean(cell2mat(P0(Iin & Iout)));
            CrossA(jj,ii) = p0_mean(1);
            CrossB(jj) = CrossB(jj)+p0_mean(2);
        end        
    end
end

if ~exist(save_folder)
    mkdir(save_folder);
end
save([save_folder,save_name,old_new{I_old_new0},mismatch_add,intercept_add,mat_tail],'CrossA','CrossB','all_color0','P0','Cin','Cout')

