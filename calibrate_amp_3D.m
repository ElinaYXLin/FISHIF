clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to calibrate between different amplifications %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNA2list.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
mask_folder = 'masks/';
mask_name = 'mask.mat';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
output_tail = '.xls';
figure_tail = '.fig';
figure_tail2 = '.png';
mat_tail = '.mat';
seg_add = '_seg';
th_add = '_th';
foci_add = '_foci';
fish_add = '_fish';
num_add = '_num';
int_add = '_int';
cmp_add = '_cmp';
sel_add = '_sel';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
RNA_add = '_RNA';
signal2_add = '_RNA2';
D3_add = '_3D';
DAPI_add = '_DAPI';
noise_add = '_noise';
mismatch_add = '_mismatch';
z_add = '_z';
xy_add= '_xy';
fit_add2 = '_fit';
FISH_add = '_FISH';
single_add = '_single';
marker_add = '_marker';
pixel_add = '_pixel';
amp_add = '_amp';
Nbin0 = 15;

channel_in = 2;
channel_ref = 1;
amp_ref = 770;

corr_add = '_corr';
cross_add = '_cross';
mismatch_add = '_0mis';
intercept_add = '_0int';

sub_dim = [2,5];
xlim0 = [0,3e3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    channel_name = eval(folder_list{list_I,5});
    [M1,M2] = size(sub_list);
    
    im_all = zeros(0);
    ref_all = zeros(0);
    amp_all = zeros(0);
    
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
%         load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'max_image');
        
        all_color = eval(folder_list{list_I,5});
        in_color = all_color{channel_in};
        k0 = strfind(sub_list{list_J,3},in_color);
        
        if ~isempty(k0)
            ind0 = k0+length(in_color)+1;
            k1 = find(sub_list{list_J,3}(ind0:end) == '_' | sub_list{list_J,3}(ind0:end) == '/' | sub_list{list_J,3}(ind0:end) == '\',1);
            ind1 = ind0+k1-2;
            amp_all = cat(1,amp_all,str2num(sub_list{list_J,3}(ind0:ind1)));
            
            fname1 = dir([image_folder,image_type]);
            temp0 = zeros(0,'uint16');
            temp1 = zeros(0,'uint16');
            for I_layer = 1:length(fname1)
                im0 = imread([image_folder,fname1(I_layer).name]);
                temp0 = cat(3,temp0,im0(:,:,channel_in)); 
                temp1 = cat(3,temp1,im0(:,:,channel_ref)); 
            end
            im_all = cat(3,im_all,mean(temp0,3));
            ref_all = cat(3,ref_all,mean(temp1,3));
        end
    end
    
    ref_mean = uint16(mean(im_all,3));
    mask0 = im2bw(ref_mean,1*graythresh(ref_mean));
    mask0 = imerode(mask0,strel('disk',10));
%     mask0 = imdilate(mask0,strel('disk',15));
    
    amp_value = amp_all;
    im_ref = im_all(:,:,find(amp_all == amp_ref,1));
% %     amp_value = unique(amp_all);
% %     im_ref = mean(im_all(:,:,amp_all == amp_ref),3);
    x0 = double(im_ref(mask0));
    
    p_value = zeros(length(amp_value),2);
    h0 = figure; maximize(h0);
    for I_amp = 1:length(amp_value)
        im1 = im_all(:,:,I_amp);
%         im1 = mean(im_all(:,:,amp_all == amp_value(I_amp)),3);
        y0 = double(im1(mask0));
        p0 = [std(y0)/std(x0),0];
%         p0 = polyfit(x0,y0,1);
        p_value(I_amp,:) = p0;

        subplot(sub_dim(1),sub_dim(2),I_amp)
            plot(x0,y0,'.','DisplayName','Data')
            hold on
            plot(xlim0,polyval(p0,xlim0),'r','DisplayName',['Fit: y = ',num2str(p0(1)),' * x + ',num2str(p0(2))])
            xlabel(['Amp_ref = ',num2str(amp_ref)])
            ylabel(['Amp = ',num2str(amp_value(I_amp))])
            title([in_color,', Amp = ',num2str(amp_value(I_amp)),', <y>/<x> = ',num2str(mean(y0)/mean(x0))])
            xlim(xlim0)
            ylim(xlim0)
            legend('show')
    end
    
    h1 = figure;
%     subplot(sub_dim(1),sub_dim(2),I_amp+1)
    subplot(1,2,1)
        plot(amp_value,p_value(:,1),'.','DisplayName','Data')
        hold on
        p_amp1 = polyfit(amp_value,p_value(:,1),2);
        plot(amp_value,polyval(p_amp1,amp_value),'r','DisplayName',['Fit: y = ',num2str(p_amp1(1)),' * x^2 + ',num2str(p_amp1(2)),' * x + ',num2str(p_amp1(3))])
        xlabel('Amp')
        ylabel('p1')
        title([in_color,', Amp_ref = ',num2str(amp_ref)])
        legend('show')
    subplot(1,2,2)
        plot(amp_value,p_value(:,2),'.','DisplayName','Data')
        hold on
        p_amp2 = polyfit(amp_value,p_value(:,2),1);
        plot(amp_value,polyval(p_amp2,amp_value),'r','DisplayName',['Fit: y = ',num2str(p_amp2(1)),' * x + ',num2str(p_amp2(2))])
        xlabel('Amp')
        ylabel('p2')
        title([in_color,', Amp_ref = ',num2str(amp_ref)])
        legend('show')
                
%%% Output: %%%============================================================
    saveas(h1,[result_folder,in_color,amp_add,intercept_add,figure_tail]);
    saveas(h1,[result_folder,in_color,amp_add,intercept_add,figure_tail2]);
    save([result_folder,in_color,amp_add,intercept_add,mat_tail],'amp_value','p_value','channel_in','in_color','p_amp1','p_amp2');
end
