clear all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_folder = 'Calibration3_FA_08112019/';
sub_folder = 'stacks/';
lsm_type = '*.lsm';
out_folder = 'Results/';
figure_tail = '.tif';
match_file = 'matchlist.xls';
standard_file = 'standard_60X.mat';
% fit_list = {'DAPI';'A488';'TMR';'A633';'A647'};
fit_list = {'DAPI';'A488';'TMR';'A647'};
list_N = length(fit_list);
record_N = zeros(1,length(fit_list));
n0 = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%% LSM file loading/resave: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lsm_temp,lsm_name] = xlsread([input_folder,sub_folder,match_file]);
    if exist([input_folder,out_folder]) ~= 7
        mkdir([input_folder,out_folder]);
    end
    
    lsm_name = cat(1,lsm_name,fit_list(1));
    N_end = find(strcmp(fit_list{1},lsm_name),1)-1;
    lsm_name = lsm_name(1:N_end);
    lsm_num0 = lsm_temp(:,1);
    output_I = 0;
    output_switch = false;

    for I_file = 1:length(lsm_name)
        input_name = lsm_name{I_file};
        input_image0 = imread([input_folder,sub_folder,input_name]); %%% Lsm file loading

        N_dim = size(input_image0);
        [X,Y] = meshgrid(1:N_dim(1),1:N_dim(2));
        XY = [reshape(X,prod(N_dim),1),reshape(Y,prod(N_dim),1)];
        Z0 = reshape(input_image0,prod(N_dim),1);
        fit2D = polyfitn(double(XY),double(Z0),n0);
        Z1 = uint16(polyvaln(fit2D,double(XY)));
        output_image0 = reshape(Z1,N_dim(1),N_dim(2));
        lsm_num(I_file,:) = [fit2D.Coefficients,fit2D.RMSE,];
%         imwrite(output_image0,[input_folder,out_folder,input_name])
        tiffwrite0(output_image0,[input_folder,out_folder,input_name])
        
        for list_I = 1:list_N
            if ~isempty(strfind(input_name,fit_list{list_I}))
                if exist(fit_list{list_I}) ~= 1
                    eval([fit_list{list_I},'=zeros(N_dim);'])
                end
                record_N(list_I) = record_N(list_I)+1;
                eval([fit_list{list_I},'=',fit_list{list_I},'+double(input_image0)/max(double(Z1));'])
                resolution(list_I) = lsm_num0(I_file);
            end
        end
    end

    for list_I = 1:list_N
        eval(['input_image0 = ',fit_list{list_I},'/',num2str(record_N(list_I)),';'])
        
        N_dim = size(input_image0);
        [X,Y] = meshgrid(1:N_dim(1),1:N_dim(2));
        XY = [reshape(X,prod(N_dim),1),reshape(Y,prod(N_dim),1)];
        Z0 = reshape(input_image0,prod(N_dim),1);
        fit2D = polyfitn(double(XY),double(Z0),n0);
        Z1 = uint16(polyvaln(fit2D,double(XY))*65535);
        output_image0 = reshape(Z1,N_dim(1),N_dim(2));
        lsm_num(length(lsm_name)+list_I,:) = [fit2D.Coefficients,fit2D.RMSE];
        lsm_num0(length(lsm_name)+list_I,:) = resolution(list_I);
%         imwrite(output_image0,[input_folder,out_folder,fit_list{list_I},figure_tail])
        tiffwrite0(output_image0,[input_folder,out_folder,fit_list{list_I},figure_tail])
        eval([fit_list{list_I},'_parameter.fit = fit2D;'])
        eval([fit_list{list_I},'_parameter.resolution = resolution(list_I);'])
        eval([fit_list{list_I},'_parameter.lattice = N_dim;'])
        eval([fit_list{list_I},'_parameter.image0 = input_image0;'])
        if exist([input_folder,out_folder,standard_file]) ~= 2
            save([input_folder,out_folder,standard_file],[fit_list{list_I},'_parameter']);
        else
            save([input_folder,out_folder,standard_file],[fit_list{list_I},'_parameter'],'-append');
        end
    end
    xlswrite([input_folder,sub_folder,match_file],cat(2,cat(1,lsm_name,fit_list),num2cell(lsm_num0),num2cell(lsm_num)));

