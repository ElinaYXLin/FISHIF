function I_true = check_embryo(embryo_name,path0)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to check whether an embryo is good to use or not %%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = fullfile('stacks/');
input_name = 'matchlist.xls';

I_true = false(size(embryo_name));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% check whether the embryo is good or not: %%%%%%%%%%%%%%%%%%%%%%
for ii = 1:length(embryo_name)
    embryo_name0 = fullfile(embryo_name{ii});
    if ~any(embryo_name0 == ':')
        embryo_name0 = [path0,embryo_name0];
    end
    I0 = strfind(embryo_name0,in_folder);
    data_folder = embryo_name0(1:I0-1);
    embryo_folder = embryo_name0(I0+length(in_folder):end);
    
    [~,~,sub_all] = xlsread([data_folder,in_folder,input_name]);
    I_file = strcmp(cellfun(@fullfile,sub_all(:,3),'UniformOutput',false),fullfile(embryo_folder));
    
    if size(sub_all,2) >= 17 && sub_all{I_file,end} == 'F'
        I_true(ii) = false;
    else
        I_true(ii) = true;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
