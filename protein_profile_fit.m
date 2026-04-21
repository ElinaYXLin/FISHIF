input_folder = '12152011/';
sub_folder = 'Results/';
input_mark = '*_protein_nucleus';
input_type = '.xls';
output_mark = '_refine';
fit_name = 'fit_parameter';

file_name = dir([input_folder,sub_folder,input_mark,output_mark,input_type]);
fit_parameter = zeros(length(file_name),3);
for I_file = 1:length(file_name)
    fname = file_name(I_file).name;
    data = xlsread([input_folder,sub_folder,fname]);
    [fit_x,fit_y,fit_parameter(I_file,:)] = protein_refine(data,[input_folder,sub_folder,fname]);
    xlswrite([input_folder,sub_folder,fname],[fit_x,fit_y]);    
end
xlswrite([input_folder,fit_name,input_type],fit_parameter);    