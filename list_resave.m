clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_folder = '01092011/';
lsm_type = '*.lsm';
out_folder = 'stacks/';
tif_name = 'stack';
figure_tail = '.tif';
match_file = 'matchlist.xls';

WGA_channel = 3;
DAPI_channel = 1;
signal1_channel = 2;
signal2_channel = 2;
Nbin = 2;
Mdim = 2;

%else_name = '10312010_3_FISH_IF_001.lsm'; 
%eWGA_channel = 4;
%eDAPI_channel = 1;
%esignal_channel = 2;
%eNbin = 2;
%eMdim = 2;

match_channel = WGA_channel;
compare_ratio = 1;
match_list = [3,3,2,4,5,3,5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% LSM file loading/resave: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsm_name = dir([input_folder,lsm_type]);
if exist([input_folder,out_folder]) ~= 7
    mkdir([input_folder,out_folder]);
end

file_name = cell(length(lsm_name)/2,2);
file_num = zeros(length(lsm_name)/2,10);
output_I = 0;
output_switch = false;

for I_file = 1:length(lsm_name)
    input_name = lsm_name(I_file).name;
    output_name = [input_name(1:(find(input_name == '.',1,'last')-1)),'/'];

    
%%% matchlist.xls output: %%%==============================================
        output_I = output_I+(~output_switch);
        file_name{output_I,output_switch+1} = output_name;
        %match_layer = ceil(length(lsm_stack)/2);
        match_layer = match_list(output_I);
        resolution = 0.091568;
        file_num(output_I,:) = [match_layer,Nbin,Mdim,match_channel,compare_ratio,WGA_channel,DAPI_channel,signal1_channel,resolution,signal2_channel];
        if ismember(output_name(length(output_name)-1),'AaBb')
            output_switch = ~output_switch;
        end
%%% =======================================================================
end

%file_num(strmatch(else_name,{lsm_name.name}),:) = [match_layer,eNbin,eMdim,eWGA_channel,compare_ratio,eWGA_channel,eDAPI_channel,esignal_channel,resolution];
xlswrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
