function [enrich_ex0,enrich_ex,legend_text] = Dualprocess5_total(total_name,iir,N_bin)

% clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real_name = 'hb';
control_name = 'nullo';
total_folder = 'enrichment_result/';
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_nullo/';
fit_folder2 = 'Histogram_default/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
output_tail = '.xls';
figure_tail = '.fig';
mat_tail = '.mat';
seg_add = '_seg';
th_add = '_th';
foci_add = '_foci';
fish_add = '_fish';
num_add = '_num';
int_add = '_int';
cmp_add = '_cmp';
sel_add = '_sel';
protein_add = '_protein';
RNA_add = '_RNA';
signal2_add = ['_',control_name];
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
reg_add = '_regulation';
fluc_add = '_fluc';
bino_add = '_bino';
local_add = '_local';
hist_add = '_hist';
fake_add = '_fake';
abs_add = '_abs';
rad_add = '_rad';
EL_add = '_EL';
D3_add = '_3D';
compare_add = '_compare';
DAPI_add = '_DAPI';
mean_add = '_mean';
slope_add = '_slope';
average_add = '_average';
spot_add = '_spot';
TX_add = '_TX';
var_add = '_var';
% sub_pos = [3,3];
sub_pos = [5,9];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [9,10,11,12,13];
enrich_out = cell(0);
enrich_out2 = cell(0);
enrich_out0 = cell(0);
enrich_out20 = cell(0);
var_out = cell(0);
var_out0 = cell(0);
var_out2 = cell(0);
var_out20 = cell(0);
TX_out = cell(0);
TX_out0 = cell(0);
TX_out2 = cell(0);
TX_out20 = cell(0);
eTX_out = cell(0);
eTX_out0 = cell(0);
eTX_out2 = cell(0);
eTX_out20 = cell(0);
EL_out = cell(0);
EL_out0 = cell(0);
EL_out2 = cell(0);
EL_out20 = cell(0);
N_out = zeros(0);
NC_out = zeros(0);
data_name = cell(0);
N_thresh0 = 50;   %%% foci number threshold for enrichment averaging
cycle_choose = [9,10,11,12,13];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('ba',folder_list(:,6)) | strcmpi('b',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3h',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
figure(111);clf
figure(576);clf
figure(5761);clf
figure(581);clf
figure(582);clf
figure(1011);clf
figure(5076);clf
figure(50761);clf
figure(50762);clf
figure(5077);clf
figure(5078);clf
figure(5081);clf
figure(5082);clf
figure(101);clf
figure(102);clf
figure(103);clf
figure(104);clf
figure(105);clf
figure(106);clf
figure(107);clf
figure(108);clf
figure(109);clf
figure(110);clf

maximize(111)
maximize(576)
maximize(5761)
maximize(581)
maximize(582)
maximize(1011)
maximize(5076)
maximize(50761)
maximize(50762)
maximize(5077)
maximize(5078)
maximize(5081)
maximize(5082)
maximize(101)
maximize(102)
maximize(103)
maximize(104)
maximize(105)
maximize(106)
maximize(107)
maximize(108)
maximize(109)
maximize(110)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN = 0;
folder_all = cell(0);

%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list, sub_all] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    if size(sub_all,2) >= 17
        I_se = ~strcmp('F',sub_all(:,17));
%         I_se = strcmp('T',sub_all(:,17));
        sub_num = sub_num(I_se,:);
        sub_list = sub_list(I_se,:);
    end
    [M1,M2] = size(sub_list);

    for list_J = 1:M1
        NN = NN+1;
        if isempty(strfind(sub_list{list_J,3},'_60X'))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
        else
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
        end
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        folder_all = cat(1,folder_all,{image_folder});
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,12);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
        clear mask_stack
        
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
%         load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
%         Inten_thresh2 = b*N_thresh;   %%% set foci intensity threshold
%         N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
%         N_cycle = round(log2(max(mask_stack(:))))+3;   %%% Calculate the nuclear cycle number
        N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};
        signal2_color = all_color{signal2_channel};
        protein_signal2_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};

        p = polyfit(nucleus_protein_profile(:,end-1),sqrt(nucleus_protein_profile(:,end)),1);
        k_DAPI = p(1);
        
        figure(111)
            imshow(double(max_image(:,:,DAPI_channel))/double(max(max(max_image(:,:,DAPI_channel)))))
            title([image_folder,': DAPI image'],'Interpreter','none')
        

%        [signal_stack,RNA_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch);   %%% load 3D image stacks
%        [nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p] = protein_profile3D(max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        I_cycle = find(N_cycle == cycle_range);
        if N_cycle > max(cycle_range)
            I_cycle = find(max(cycle_range) == cycle_range);
        elseif N_cycle < min(cycle_range)
            I_cycle = find(min(cycle_range) == cycle_range);
        end
        if ~isempty(I_cycle)
            cycle_pos0(I_cycle) = cycle_pos0(I_cycle)+1;
            sub_pos(3) = sub2ind(sub_pos([2,1]),cycle_pos0(I_cycle),I_cycle);
        else
            sub_pos(3) = 0;
        end
        
        figure(1011)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            imshow(double(max_image(:,:,DAPI_channel))/double(max(max(max_image(:,:,DAPI_channel)))))
            title([image_folder,char(10),'DAPI image, cycle = ',num2str(N_cycle),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI)],'Interpreter','none')
       
        
%         [enrich_temp,enrich_temp2,enrich_temp0,enrich_temp20,var_temp,var_temp0,N_temp] = enrichment_compare(nucleus_protein_profile,nucleus_RNA_profile/b,foci_data,fake_data,r_size,foci_data2,fake_data2,signal2_add(2:end),image_folder,N_cycle,sub_pos,k_DAPI,true,iir,z_size,real_name);
        [enrich_temp,enrich_temp2,enrich_temp0,enrich_temp20,var_temp,var_temp2,var_temp0,var_temp20,TX_temp,TX_temp2,TX_temp0,TX_temp20,eTX_temp,eTX_temp2,eTX_temp0,eTX_temp20,EL_temp,EL_temp2,EL_temp0,EL_temp20,N_temp] = enrichment_compare(nucleus_protein_profile,nucleus_RNA_profile,foci_data,fake_data,h,r_size,N_thresh,foci_data2,fake_data2,signal2_add(2:end),['Sample ',num2str(NN)],N_cycle,sub_pos,k_DAPI,N_bin,true,iir,z_size,real_name);
        enrich_out = cat(1,enrich_out,{enrich_temp});
        enrich_out2 = cat(1,enrich_out2,{enrich_temp2});
        enrich_out0 = cat(1,enrich_out0,{enrich_temp0});
        enrich_out20 = cat(1,enrich_out20,enrich_temp20);
        var_out = cat(1,var_out,{var_temp});
        var_out0 = cat(1,var_out0,{var_temp0});
        var_out2 = cat(1,var_out2,{var_temp2});
        var_out20 = cat(1,var_out20,var_temp20);
        TX_out = cat(1,TX_out,{TX_temp});
        TX_out0 = cat(1,TX_out0,{TX_temp0});
        TX_out2 = cat(1,TX_out2,{TX_temp2});
        TX_out20 = cat(1,TX_out20,TX_temp20);
        eTX_out = cat(1,eTX_out,{eTX_temp});
        eTX_out0 = cat(1,eTX_out0,{eTX_temp0});
        eTX_out2 = cat(1,eTX_out2,{eTX_temp2});
        eTX_out20 = cat(1,eTX_out20,eTX_temp20);
        EL_out = cat(1,EL_out,{EL_temp});
        EL_out0 = cat(1,EL_out0,{EL_temp0});
        EL_out2 = cat(1,EL_out2,{EL_temp2});
        EL_out20 = cat(1,EL_out20,EL_temp20);
        N_out = cat(1,N_out,N_temp);
        NC_out = cat(1,NC_out,N_cycle);
        data_name = cat(1,data_name,{image_folder});
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        
        saveas(111,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),DAPI_add,figure_tail]);
%        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
        saveas(576,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,compare_add,figure_tail]);
        saveas(5761,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,TX_add,D3_add,compare_add,figure_tail]);
        saveas(581,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,compare_add,figure_tail]);
        saveas(582,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,EL_add,D3_add,compare_add,figure_tail]);
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','foci_signal2_profile','cytoplasmic_signal2_profile','foci_data2','fake_data2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end
%     xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));



%% Data output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if total_name(end) == '/' || total_name(end) == '\'
    total_name = total_name(1:end-1);
end
result_folder2 = [total_folder,total_name,'/'];
total_name0 = total_name;
total_name0((total_name == '/') | (total_name == '\')) = '_';

if exist(result_folder2) ~= 7
    mkdir(result_folder2);
end


saveas(1011,[result_folder2,total_name0,DAPI_add,compare_add,figure_tail]);
saveas(5076,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,compare_add,figure_tail]);
saveas(50761,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,TX_add,D3_add,compare_add,figure_tail]);
saveas(50762,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,var_add,D3_add,compare_add,figure_tail]);
saveas(5077,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,mean_add,D3_add,compare_add,figure_tail]);
saveas(5078,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,slope_add,D3_add,compare_add,figure_tail]);
saveas(5081,[result_folder2,total_name0,protein_add,local_add,abs_add,rad_add,D3_add,compare_add,figure_tail]);
saveas(5082,[result_folder2,total_name0,protein_add,local_add,abs_add,EL_add,D3_add,compare_add,figure_tail]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Average over all embryos: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_bin = 25;
unit1 = '#';
unit2 = 'M';
color_code = mycolors(size(enrich_out{1},2)/2);
color_code = [color_code;color_code(1:end-1,:)];
cycle_sel = cat(2,num2cell(cycle_range),{cycle_choose});
name_sel = cat(2,cellstr([num2str(cycle_range'),repmat('/',length(cycle_range),1)])',{''});

for I_sel = 1:length(cycle_sel)
%     cycle_sel{I_sel}
    good_data = (N_out >= N_thresh0) & ismember(NC_out,cycle_sel{I_sel});
    good_data0 = ismember(NC_out,cycle_sel{I_sel});
    if ~isempty(good_data)
        enrich_out00 = [cell2mat(enrich_out(good_data)),cell2mat(enrich_out0(good_data))];
        var_out00 = [cell2mat(var_out(good_data)),cell2mat(var_out0(good_data))];
        TX_out00 = [cell2mat(TX_out(good_data)),cell2mat(TX_out0(good_data))];
        eTX_out00 = [cell2mat(eTX_out(good_data)),cell2mat(eTX_out0(good_data))];
        EL_out00 = [cell2mat(EL_out(good_data)),cell2mat(EL_out0(good_data))];
        enrich_embryo = [cell2mat(enrich_out(good_data0)),cell2mat(enrich_out0(good_data0))];
        var_embryo = [cell2mat(var_out(good_data0)),cell2mat(var_out0(good_data0))];
        TX_embryo = [cell2mat(TX_out(good_data0)),cell2mat(TX_out0(good_data0))];
        eTX_embryo = [cell2mat(eTX_out(good_data0)),cell2mat(eTX_out0(good_data0))];
        EL_embryo = [cell2mat(EL_out(good_data0)),cell2mat(EL_out0(good_data0))];
        enrich_ex0 = enrich_out00;

        legend_text = {[real_name,' - background'],['fake ',real_name,' - background'],[real_name,' - fake'],[control_name,' - background'],[control_name,' fake - background'],[control_name,' - fake'],[real_name,' - ',control_name], ...
                       [real_name,' - background all'],['fake ',real_name,' - background all'],[real_name,' - fake all'],[control_name,' - background all'],[control_name,' fake - background all'],[control_name,' - fake all']};
        legend_textv = {[real_name,' - background'],['fake ',real_name,' - background'],['(fake ',real_name,' - background) fit'],[real_name,' - (fake fit)'],[control_name,' - background'],[control_name,' fake - background'],['(',control_name,' fake - background) fit'],[control_name,' - (fake fit)'],[real_name,' - ',control_name], ...
                        [real_name,' - background all'],['fake ',real_name,' - background all'],['(fake ',real_name,' - background all) fit'],[real_name,' - (fake all fit)'],[control_name,' - background all'],[control_name,' fake - background all'],['(',control_name,' fake - background all) fit'],[control_name,' - (fake all fit)']};
        diff_index = [3,6,10,13];

        figure(101)
        clf
        for I_plot = 1:size(enrich_out00,2)/2
            enrich_out0temp = enrich_out00(:,(I_plot*2-1):(I_plot*2));
            [xtemp1,IX1] = sort(enrich_out0temp(:,1));
            ytemp1 = enrich_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end
        %     p_fit = polyfit(x1_f,y1_f,1);
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        %     plot(x1_f,p_fit(1)*x1_f+p_fit(2),[color_code(I_plot),'--'])
        end

        title([total_name,char(10),'Mean foci enrichment vs protein concentration (average over ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
        ylabel(['Protein enrichment (',unit1,')'])
        xlabel(['P_m_e_a_n (',unit2,')'])
        legend(legend_text)


        figure(103)
        clf
        enrich_mean = zeros(1,size(enrich_out00,2)/2);
        enrich_err  = zeros(1,size(enrich_out00,2)/2);
        for I_plot = 1:size(enrich_out00,2)/2
            enrich_mean(I_plot) = mean(enrich_out00(:,(I_plot*2)));
            enrich_err(I_plot)  = std0(enrich_out00(:,(I_plot*2)));
            bar(I_plot,enrich_mean(I_plot),'FaceColor',color_code(I_plot,:));
            hold on
        end
        errorbar([1:size(enrich_out00,2)/2],enrich_mean,enrich_err,'k.')

        title([total_name,char(10),'Comparison of enrichment calculation  (average over ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
        ylabel(['Protein enrichment (',unit1,')'])
        % set(gca,'XTick',[1:length(legend_text)],'XTickLabel',legend_text)
        xticklabel_rotate([1:length(legend_text)],45,legend_text,'FontSize',14,'FontWeight','bold')


        figure(105)
        clf
        for I_plot = 1:size(var_out00,2)/2
            var_out0temp = var_out00(:,(I_plot*2-1):(I_plot*2));
            [xtemp1,IX1] = sort(var_out0temp(:,1));
            ytemp1 = var_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        end

        title([total_name,char(10),'Foci enrichment variance vs protein concentration (average over ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
        ylabel(['Protein enrichment variance (',unit1,')'])
        xlabel(['P_m_e_a_n (',unit2,')'])
        legend(legend_text)



        figure(107)
        clf
        for I_plot = 1:size(TX_out00,2)/3
            TX_out0temp = TX_out00(:,(I_plot*3-2):(I_plot*3));
            [ztemp1,IX1] = sort(TX_out0temp(:,3));
            xtemp1 = TX_out0temp(IX1,1);
            ytemp1 = TX_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            z1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            z1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_f(I_bin) = mean(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_ferr(I_bin) = std0(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        end
        for I_plot = 1:size(eTX_out00,2)/3
            eTX_out0temp = eTX_out00(:,(I_plot*3-2):(I_plot*3));
            ztemp1 = eTX_out0temp(:,3);
            xtemp1 = eTX_out0temp(:,1);
            ytemp1 = eTX_out0temp(:,2);
            Igood = ~isnan(xtemp1);
            ztemp1 = ztemp1(Igood);
            xtemp1 = xtemp1(Igood);
            ytemp1 = ytemp1(Igood);
            
            EL_all = unique(ztemp1);
            x1_f = zeros(size(EL_all));
            y1_f = zeros(size(EL_all));
            z1_f = EL_all;
            x1_ferr = zeros(size(EL_all));
            y1_ferr = zeros(size(EL_all));
            z1_ferr = zeros(size(EL_all));
            for I_bin = 1:length(EL_all)
                x1_f(I_bin) = mean(xtemp1(ztemp1 == EL_all(I_bin)));
                x1_ferr(I_bin) = std0(xtemp1(ztemp1 == EL_all(I_bin)));
                y1_f(I_bin) = mean(ytemp1(ztemp1 == EL_all(I_bin)));
                y1_ferr(I_bin) = std0(ytemp1(ztemp1 == EL_all(I_bin)));
            end
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        end

        title([total_name,char(10),'Mean TX level vs mean foci enrichment (average over ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
        ylabel(['Mean TX level (',unit1,')'])
        xlabel(['Mean foci enrichment (',unit1,')'])
        legend_text_107 = {['(',real_name,' - background) bin'],['(',real_name,' - fake) bin'],['(',control_name,' - background) bin'],['(',control_name,' - fake) bin'],['(',real_name,' - ',control_name,') bin'],['(',real_name,' - background) all bin'],['(',real_name,' - fake) all bin'],['(',control_name,' - background) all bin'],['(',control_name,' - fake) all bin']};
        legend_text_107 = cat(2,legend_text_107,{['(',real_name,' - background) bin (EL)'],['(',real_name,' - fake) bin (EL)'],['(',control_name,' - background) bin (EL)'],['(',control_name,' - fake) bin (EL)'],['(',real_name,' - ',control_name,') bin (EL)'],['(',real_name,' - background) all bin (EL)'],['(',real_name,' - fake) all bin (EL)'],['(',control_name,' - background) all bin (EL)'],['(',control_name,' - fake) all bin (EL)']});
        legend(legend_text_107)


        figure(109)
        clf
        for I_plot = 1:size(EL_out00,2)/2
            xtemp1 = EL_out00(:,(I_plot*2-1));
            ytemp1 = EL_out00(:,(I_plot*2));
            x1_f = unique(xtemp1);
            y1_f = zeros(size(x1_f));
            y1_ferr = zeros(size(x1_f));
            for I_bin = 1:length(x1_f)
                I_temp = (xtemp1 == x1_f(I_bin)) & ~isnan(ytemp1);
                y1_f(I_bin) = mean(ytemp1(I_temp));
                y1_ferr(I_bin) = std0(ytemp1(I_temp));
            end
            errorbar(x1_f,y1_f,y1_ferr,'Color',color_code(I_plot,:),'LineWidth',2)
            hold on
        end

        title([total_name,char(10),'Foci enrichment vs EL (average over ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
        ylabel(['Protein enrichment (',unit1,')'])
        xlabel(['EL'])
        legend(legend_text)

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        %% Average over all foci: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        enrich_out00 = cell2mat(enrich_out2(good_data));
        enrich_out001 = cell2mat(enrich_out20(good_data,1));
        enrich_out002 = cell2mat(enrich_out20(good_data,2));
        var_out00 = cell2mat(var_out2(good_data));
        var_out001 = cell2mat(var_out20(good_data,1));
        var_out002 = cell2mat(var_out20(good_data,2));
        TX_out00 = cell2mat(TX_out2(good_data));
        TX_out001 = cell2mat(TX_out20(good_data,1));
        TX_out002 = cell2mat(TX_out20(good_data,2));
        eTX_out00 = cell2mat(eTX_out2(good_data));
        eTX_out001 = cell2mat(eTX_out20(good_data,1));
        eTX_out002 = cell2mat(eTX_out20(good_data,2));
        EL_out00 = cell2mat(EL_out2(good_data));
        EL_out001 = cell2mat(EL_out20(good_data,1));
        EL_out002 = cell2mat(EL_out20(good_data,2));
        enrich_foci = {cell2mat(enrich_out2(good_data0)),cell2mat(enrich_out20(good_data0,1)),cell2mat(enrich_out20(good_data0,2))};
        var_foci = {cell2mat(var_out2(good_data0)),cell2mat(var_out20(good_data0,1)),cell2mat(var_out20(good_data0,2))};
        TX_foci = {cell2mat(TX_out2(good_data0)),cell2mat(TX_out20(good_data0,1)),cell2mat(TX_out20(good_data0,2))};
        eTX_foci = {cell2mat(eTX_out2(good_data0)),cell2mat(eTX_out20(good_data0,1)),cell2mat(eTX_out20(good_data0,2))};
        EL_foci = {cell2mat(EL_out2(good_data0)),cell2mat(EL_out20(good_data0,1)),cell2mat(EL_out20(good_data0,2))};
        
        enrich_out00 = enrich_out00(:,1:end-1);
        enrich_out001 = enrich_out001(:,1:end-1);
        enrich_out002 = enrich_out002(:,1:end-1);
        var_out00 = var_out00(:,1:end-1);
        var_out001 = var_out001(:,1:end-1);
        var_out002 = var_out002(:,1:end-1);
        enrich_ex = {enrich_out00,enrich_out001,enrich_out002};

        figure(102)
        clf
        for I_plot = 1:size(enrich_out00,2)/2
            enrich_out0temp = enrich_out00(:,(I_plot*2-1):(I_plot*2));
            [xtemp1,IX1] = sort(enrich_out0temp(:,1));
            ytemp1 = enrich_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end
        %     p_fit = polyfit(x1_f,y1_f,1);
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        %     plot(x1_f,p_fit(1)*x1_f+p_fit(2),[color_code(I_plot),'--'])
        end
        N0 = size(enrich_out00,2)/2;
        for I_plot = 1:size(enrich_out001,2)/2
            enrich_out0temp = enrich_out001(:,(I_plot*2-1):(I_plot*2));
            [xtemp1,IX1] = sort(enrich_out0temp(:,1));
            ytemp1 = enrich_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end
        %     p_fit = polyfit(x1_f,y1_f,1);
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
            hold on
        %     plot(x1_f,p_fit(1)*x1_f+p_fit(2),[color_code(I_plot),'--'])
        end
        N0 = size(enrich_out00,2)/2+size(enrich_out001,2)/2;
        for I_plot = 1:size(enrich_out002,2)/2
            enrich_out0temp = enrich_out002(:,(I_plot*2-1):(I_plot*2));
            [xtemp1,IX1] = sort(enrich_out0temp(:,1));
            ytemp1 = enrich_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end
        %     p_fit = polyfit(x1_f,y1_f,1);
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
            hold on
        %     plot(x1_f,p_fit(1)*x1_f+p_fit(2),[color_code(I_plot),'--'])
        end

        title([total_name,char(10),'Mean foci enrichment vs protein concentration (average over ',num2str(size(enrich_out00,1)),' foci)'],'Interpreter','none')
        ylabel(['Protein enrichment (',unit1,')'])
        xlabel(['P_m_e_a_n (',unit2,')'])
        legend(legend_text)
        N0 = size(enrich_out00,2)/2+size(enrich_out001,2)/2+size(enrich_out002,2)/2;


        figure(104)
        clf
        enrich_mean = zeros(1,N0);
        enrich_err  = zeros(1,N0);
        enrich_mean = [mean(enrich_out00(:,2:2:end)),mean(enrich_out001(:,2:2:end)),mean(enrich_out002(:,2:2:end))];
        enrich_err  = [std0(enrich_out00(:,2:2:end)),std0(enrich_out001(:,2:2:end)),std0(enrich_out002(:,2:2:end))];
        if N0
            for I_plot = 1:N0
                bar(I_plot,enrich_mean(I_plot),'FaceColor',color_code(I_plot,:));
                hold on
            end
            errorbar([1:N0],enrich_mean,enrich_err,'k.')
        end

        title([total_name,char(10),'Comparison of enrichment calculation  (average over ',num2str(size(enrich_out00,1)),' foci)'],'Interpreter','none')
        ylabel(['Protein enrichment (',unit1,')'])
        % set(gca,'XTick',[1:length(legend_text)],'XTickLabel',legend_text)
        xticklabel_rotate([1:length(legend_text)],45,legend_text,'FontSize',14,'FontWeight','bold')



        figure(106)
        clf
        I_switch = [3,6,7];
        for I_plot = 1:size(var_out00,2)/3
            var_out0temp = var_out00(:,(I_plot*3-2):(I_plot*3));
            [xtemp1,IX1] = sort(var_out0temp(:,1));
            ytemp1 = var_out0temp(IX1,2);
            ytemp0 = var_out0temp(IX1,3);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                if any(I_plot == I_switch)
                    y1_f(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp0(((I_bin-1)*N_temp+1):I_bin*N_temp));
                    y1_ferr(I_bin) = y1_f(I_bin)/sqrt((N_temp)/2);
                else
                    y1_f(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp0(((I_bin-1)*N_temp+1):I_bin*N_temp));
                    y1_ferr(I_bin) = y1_f(I_bin)/sqrt((N_temp)/2);
                end
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        end

        N0 = size(var_out00,2)/3;
        for I_plot = 1:size(var_out001,2)/3
            var_out0temp = var_out001(:,(I_plot*3-2):(I_plot*3));
            [xtemp1,IX1] = sort(var_out0temp(:,1));
            ytemp1 = var_out0temp(IX1,2);
            ytemp0 = var_out0temp(IX1,3);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                if any(I_plot == I_switch)
                    y1_f(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp0(((I_bin-1)*N_temp+1):I_bin*N_temp));
                    y1_ferr(I_bin) = y1_f(I_bin)/sqrt((N_temp)/2);
                else
                    y1_f(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp0(((I_bin-1)*N_temp+1):I_bin*N_temp));
                    y1_ferr(I_bin) = y1_f(I_bin)/sqrt((N_temp)/2);
                end
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
            hold on
        end

        N0 = size(var_out00,2)/3+size(var_out001,2)/3;
        for I_plot = 1:size(var_out002,2)/3
            var_out0temp = var_out002(:,(I_plot*3-2):(I_plot*3));
            [xtemp1,IX1] = sort(var_out0temp(:,1));
            ytemp1 = var_out0temp(IX1,2);
            ytemp0 = var_out0temp(IX1,3);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                if any(I_plot == I_switch)
                    y1_f(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp0(((I_bin-1)*N_temp+1):I_bin*N_temp));
                    y1_ferr(I_bin) = y1_f(I_bin)/sqrt((N_temp)/2);
                else
                    y1_f(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp0(((I_bin-1)*N_temp+1):I_bin*N_temp));
                    y1_ferr(I_bin) = y1_f(I_bin)/sqrt((N_temp)/2);
                end
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
            hold on
        end

        title([total_name,char(10),'Foci enrichment variance vs protein concentration (average over ',num2str(size(enrich_out00,1)),' foci)'],'Interpreter','none')
        ylabel(['Protein enrichment variance (',unit1,')'])
        xlabel(['P_m_e_a_n (',unit2,')'])
        legend(legend_text)



        figure(108)
        clf
        for I_plot = 1:size(TX_out00,2)/3
            TX_out0temp = TX_out00(:,(I_plot*3-2):(I_plot*3));
            [ztemp1,IX1] = sort(TX_out0temp(:,3));
            xtemp1 = TX_out0temp(IX1,1);
            ytemp1 = TX_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            z1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            z1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_f(I_bin) = mean(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_ferr(I_bin) = std0(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        end

        N0 = size(var_out00,2)/3;
        for I_plot = 1:size(TX_out001,2)/3
            TX_out0temp = TX_out001(:,(I_plot*3-2):(I_plot*3));
            [ztemp1,IX1] = sort(TX_out0temp(:,3));
            xtemp1 = TX_out0temp(IX1,1);
            ytemp1 = TX_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            z1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            z1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_f(I_bin) = mean(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_ferr(I_bin) = std0(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
        end

        N0 = size(var_out00,2)/3+size(var_out001,2)/3;
        for I_plot = 1:size(TX_out002,2)/3
            TX_out0temp = TX_out002(:,(I_plot*3-2):(I_plot*3));
            [ztemp1,IX1] = sort(TX_out0temp(:,3));
            xtemp1 = TX_out0temp(IX1,1);
            ytemp1 = TX_out0temp(IX1,2);
            N_temp = floor(length(xtemp1)/N_bin);
            x1_f = zeros(1,N_bin);
            y1_f = zeros(1,N_bin);
            z1_f = zeros(1,N_bin);
            x1_ferr = zeros(1,N_bin);
            y1_ferr = zeros(1,N_bin);
            z1_ferr = zeros(1,N_bin);
            for I_bin = 1:N_bin
                x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_f(I_bin) = mean(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
                z1_ferr(I_bin) = std0(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for I_plot = 1:size(eTX_out00,2)/3
            eTX_out0temp = eTX_out00(:,(I_plot*3-2):(I_plot*3));
            ztemp1 = eTX_out0temp(:,3);
            xtemp1 = eTX_out0temp(:,1);
            ytemp1 = eTX_out0temp(:,2);
            Lbin = EL_all(2)-EL_all(1);
            x1_f = zeros(size(EL_all));
            y1_f = zeros(size(EL_all));
            z1_f = EL_all;
            x1_ferr = zeros(size(EL_all));
            y1_ferr = zeros(size(EL_all));
            z1_ferr = zeros(size(EL_all));
            for I_bin = 1:length(EL_all)
                x1_f(I_bin) = mean(xtemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                x1_ferr(I_bin) = std0(xtemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                y1_f(I_bin) = mean(ytemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                y1_ferr(I_bin) = std0(ytemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2)
            hold on
        end

        N0 = size(var_out00,2)/3;
        for I_plot = 1:size(eTX_out001,2)/3
            eTX_out0temp = eTX_out001(:,(I_plot*3-2):(I_plot*3));
            ztemp1 = eTX_out0temp(:,3);
            xtemp1 = eTX_out0temp(:,1);
            ytemp1 = eTX_out0temp(:,2);
            Lbin = EL_all(2)-EL_all(1);
            x1_f = zeros(size(EL_all));
            y1_f = zeros(size(EL_all));
            z1_f = EL_all;
            x1_ferr = zeros(size(EL_all));
            y1_ferr = zeros(size(EL_all));
            z1_ferr = zeros(size(EL_all));
            for I_bin = 1:length(EL_all)
                x1_f(I_bin) = mean(xtemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                x1_ferr(I_bin) = std0(xtemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                y1_f(I_bin) = mean(ytemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                y1_ferr(I_bin) = std0(ytemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
        end

        N0 = size(var_out00,2)/3+size(var_out001,2)/3;
        for I_plot = 1:size(eTX_out002,2)/3
            eTX_out0temp = eTX_out002(:,(I_plot*3-2):(I_plot*3));
            ztemp1 = eTX_out0temp(:,3);
            xtemp1 = eTX_out0temp(:,1);
            ytemp1 = eTX_out0temp(:,2);
            Lbin = EL_all(2)-EL_all(1);
            x1_f = zeros(size(EL_all));
            y1_f = zeros(size(EL_all));
            z1_f = EL_all;
            x1_ferr = zeros(size(EL_all));
            y1_ferr = zeros(size(EL_all));
            z1_ferr = zeros(size(EL_all));
            for I_bin = 1:length(EL_all)
                x1_f(I_bin) = mean(xtemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                x1_ferr(I_bin) = std0(xtemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                y1_f(I_bin) = mean(ytemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
                y1_ferr(I_bin) = std0(ytemp1((ztemp1 >= EL_all(I_bin)-Lbin) & (ztemp1 <= EL_all(I_bin)+Lbin)));
            end

            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_code(I_plot+N0,:),color_code(I_plot+N0,:),'o','LineWidth',2)
        end

        title([total_name,char(10),'Mean TX level vs foci enrichment (average over ',num2str(size(enrich_out00,1)),' foci)'],'Interpreter','none')
        ylabel(['Mean TX level (',unit1,')'])
        xlabel(['Mean foci enrichment (',unit1,')'])
        legend(legend_text_107)




        figure(110)
        clf
        Lmin = 0.1;%0.15;
        Lmax = 0.75;%0.85;
        Lbin = 0.05;
        EL_all = Lmin:Lbin:Lmax;

        for I_plot = 1:size(EL_out00,2)/2
            xtemp = EL_out00(:,(I_plot*2-1));
            ytemp = EL_out00(:,I_plot*2);
            y1_f = zeros(size(EL_all));
            y1_ferr = zeros(size(EL_all));
            for I_bin = 1:length(EL_all)
                I_temp = (xtemp >= EL_all(I_bin)-Lbin) & (xtemp <= EL_all(I_bin)+Lbin);
                y1_f(I_bin) = mean(ytemp(I_temp));
                y1_ferr(I_bin) = std0(ytemp(I_temp));
            end
            errorbar(EL_all,y1_f,y1_ferr,'Color',color_code(I_plot,:),'LineWidth',2)
            hold on
        end
        N0 = size(EL_out00,2)/2;

        for I_plot = 1:size(EL_out001,2)/2
            xtemp = EL_out001(:,(I_plot*2-1));
            ytemp = EL_out001(:,I_plot*2);
            y1_f = zeros(size(EL_all));
            y1_ferr = zeros(size(EL_all));
            for I_bin = 1:length(EL_all)
                I_temp = (xtemp >= EL_all(I_bin)-Lbin) & (xtemp <= EL_all(I_bin)+Lbin);
                y1_f(I_bin) = mean(ytemp(I_temp));
                y1_ferr(I_bin) = std0(ytemp(I_temp));
            end
            errorbar(EL_all,y1_f,y1_ferr,'Color',color_code(I_plot+N0,:),'LineWidth',2)
            hold on
        end
        N0 = size(EL_out00,2)/2+size(EL_out001,2)/2;

        for I_plot = 1:size(EL_out002,2)/2
            xtemp = EL_out002(:,(I_plot*2-1));
            ytemp = EL_out002(:,I_plot*2);
            y1_f = zeros(size(EL_all));
            y1_ferr = zeros(size(EL_all));
            for I_bin = 1:length(EL_all)
                I_temp = (xtemp >= EL_all(I_bin)-Lbin) & (xtemp <= EL_all(I_bin)+Lbin);
                y1_f(I_bin) = mean(ytemp(I_temp));
                y1_ferr(I_bin) = std0(ytemp(I_temp));
            end
            errorbar(EL_all,y1_f,y1_ferr,'Color',color_code(I_plot+N0,:),'LineWidth',2)
            hold on
        end

        title([total_name,char(10),'Mean foci enrichment vs EL (average over ',num2str(size(enrich_out00,1)),' foci)'],'Interpreter','none')
        ylabel(['Protein enrichment (',unit1,')'])
        xlabel('EL')
        legend(legend_text)

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    %     if total_name(end) == '/' || total_name(end) == '\' || total_name(end) == '_'
    %         total_name = total_name(1:end-1);
    %     end
        result_folder2 = [total_folder,total_name,'/',name_sel{I_sel}];
    %     total_name((total_name == '/') | (total_name == '\')) = '_';

        if exist(result_folder2) ~= 7
            mkdir(result_folder2);
        end


        saveas(101,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,figure_tail]);
        saveas(102,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,spot_add,figure_tail]);
        saveas(103,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,mean_add,figure_tail]);
        saveas(104,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,mean_add,spot_add,figure_tail]);
        saveas(105,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,var_add,figure_tail]);
        saveas(106,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,var_add,spot_add,figure_tail]);
        saveas(107,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,TX_add,figure_tail]);
        saveas(108,[result_folder2,total_name0,protein_add,local_add,abs_add,foci_add,D3_add,average_add,TX_add,spot_add,figure_tail]);
        saveas(109,[result_folder2,total_name0,protein_add,local_add,abs_add,EL_add,D3_add,average_add,figure_tail]);
        saveas(110,[result_folder2,total_name0,protein_add,local_add,abs_add,EL_add,D3_add,average_add,spot_add,figure_tail]);
        save([result_folder2,total_name0,mat_tail],'enrich_embryo','var_embryo','TX_embryo','eTX_embryo','EL_embryo','enrich_foci','var_foci','TX_foci','eTX_foci','EL_foci','folder_all');
        try
        xlswrite([result_folder2,total_name0,input_name],cat(2,data_name,num2cell(good_data)));
        catch

        end
    end
end    
toc