function crosstalk_compare
clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_folder = 'crosstalk/';
lsm_type = '*.lsm';
out_folder = 'results/';
tif_name = 'stack';
figure_tail = '.tif';
figure_tail2 = '.fig';
fit_out = 'fitting_result';
xls_tail = '.xls';
match_file = 'matchlist.xls';
N_tiff = 20;
w_tiff = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%% LSM file loading/resave: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsm_name = dir([input_folder,lsm_type]);
if exist([input_folder,out_folder]) ~= 7
    mkdir([input_folder,out_folder]);
end

file_name = cell(length(lsm_name)/2,2);
output_I = 0;
output_switch = false;
p0 = zeros(length(lsm_name),2);

for I_file = 1:length(lsm_name)
    input_name = lsm_name(I_file).name;
    lsm_stack = tiffread([input_folder,input_name]); %%% Lsm file loading

    tiff_image = lsm_stack.data;
    tiff_A488 = tiff_image{1};
    tiff_TMR = tiff_image{2};
    tiff_A488 = double(tiff_A488(:));
    tiff_TMR = double(tiff_TMR(:));

    figure(1)
    clf
    plot(tiff_TMR,tiff_A488,'b.')
    hold on

    tbin = (max(tiff_TMR)-min(tiff_TMR))/N_tiff;
    tiff_window = (max(tiff_TMR)-min(tiff_TMR))*w_tiff;
    tiff_bin = [(min(tiff_TMR)+tbin/2):tbin:(max(tiff_TMR)-tbin/2)];
    tiff_mean = zeros(size(tiff_bin));
    tiff_err = zeros(size(tiff_bin));
    for I_bin = 1:length(tiff_bin)
        tiff_mean(I_bin) = mean(tiff_A488((tiff_TMR >= (tiff_bin(I_bin)-tiff_window))&(tiff_TMR <= (tiff_bin(I_bin)+tiff_window))));
        tiff_err(I_bin) = std0(tiff_A488((tiff_TMR >= (tiff_bin(I_bin)-tiff_window))&(tiff_TMR <= (tiff_bin(I_bin)+tiff_window))));
    end
    errorbar(tiff_bin,tiff_mean,tiff_err,'r')

    %p = polyfit(tiff_TMR,tiff_A488,1);
    p = polyfit(tiff_bin(2:12),tiff_mean(2:12),1);
    x = [min(tiff_TMR),max(tiff_TMR)];
    plot(x,x*p(1)+p(2),'g')

    xlabel('TMR channel')
    ylabel('A488 channel')
    title(['Channel crosstalk plot: ',input_name,', I_A488 = ',num2str(p(1)),' * I_TMR + ',num2str(p(2))])
    saveas(1,[input_folder,out_folder,input_name(1:(find(input_name == '.',1,'last')-1)),figure_tail2])
    p0(I_file,:) = p;
end
xlswrite([input_folder,out_folder,fit_out,xls_tail],p0)
    
    
function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);