%% Subfunction to load variable from mat file: %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out = var_load(file_name,var_name)
S = load(file_name,var_name);
var_out = eval(['S.',var_name]);


