clear all
%tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Enrichmentlist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
output_folder = 'enrichment_result/';
hist_folder = 'Histogram/';
hist_tail = '_raw.xls';
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
mark_code = '+o*.xsd^v><ph';
ccode_foci = {[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[0,0,0]};
ccode_fake = {[1,0,0]/2,[0,1,0]/2,[0,0,1]/2,[1,1,0]/2,[1,0,1]/2,[0,1,1]/2,[1,1,1]/2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, group_list] = xlsread(list_name);
hist_text = cell(0);
hist_I = 0;
hist_legend = cell(0);

for list_G = 1:size(group_list,1)
    group_name = group_list{list_G,1};
    folder_list = eval(group_list{list_G,2});
    %folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
    N1 = length(folder_list);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data73a = zeros(0);
    data73b = zeros(0);
    data74  = zeros(0);
    data75  = zeros(0);
    data76a = zeros(0);
    data76b = zeros(0);
    data77  = zeros(0);
    data_legend0 = cell(0);
    data_legend  = cell(0);
    data_legend1 = cell(0);
    
    h = cell(1,10);
    foci_data0 = cell(size(h));
    fake_data0 = cell(size(h));
    p_ratio00 = cell(size(h));
    fp_ratio00 = cell(size(h));
    op_ratio00 = cell(size(h));
    ofp_ratio00 = cell(size(h));
    for ir = 1:length(h)
        foci_data0{ir} = zeros(0);
        fake_data0{ir} = zeros(0);
        p_ratio00{ir} = zeros(0);
        fp_ratio00{ir} = zeros(0);
        op_ratio00{ir} = zeros(0);
        ofp_ratio00{ir} = zeros(0);
    end

    mean_p   = zeros(0);
    mean_fp  = zeros(0);
    mean_op  = zeros(0);
    mean_ofp = zeros(0);
    std_p   = zeros(0);
    std_fp  = zeros(0);
    std_op  = zeros(0);
    std_ofp = zeros(0);

    clear t_absolute
    %% Enrichment analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for list_I = 1:N1
        [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
        [M1,M2] = size(sub_list);
        I_mark = 0;
        
        for list_J = 1:M1
            I_mark = I_mark+1;
            image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
            result_folder = [folder_list{list_I,1},out_folder];

            load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
            if ~exist('t_absolute')
                t_absolute = false;
            end
            [out73a,out73b,out74,out75,out76a,out76b,out77,p_ratio0,fp_ratio0,op_ratio0,ofp_ratio0,r_size,legend_text0,legend_text,legend_text1] = enrichment_profile(foci_data,fake_data,h,r_size,image_folder,mark_code(I_mark));
            temp_p_mean = zeros(size(h));
            temp_fp_mean = zeros(size(h));
            temp_op_mean = zeros(size(h));
            temp_ofp_mean = zeros(size(h));
            temp_p_std = zeros(size(h));
            temp_fp_std = zeros(size(h));
            temp_op_std = zeros(size(h));
            temp_ofp_std = zeros(size(h));

            for ir = 1:length(h)
                foci_data0{ir} = cat(1,foci_data0{ir},foci_data{ir});
                fake_data0{ir} = cat(1,fake_data0{ir},fake_data{ir});
                p_ratio00{ir}   = cat(1,p_ratio00{ir},p_ratio0{ir});
                fp_ratio00{ir}  = cat(1,fp_ratio00{ir},fp_ratio0{ir});
                op_ratio00{ir}  = cat(1,op_ratio00{ir},op_ratio0{ir});
                ofp_ratio00{ir} = cat(1,ofp_ratio00{ir},ofp_ratio0{ir});
                
                temp_p_mean(ir) = mean(p_ratio0{ir});
                temp_fp_mean(ir) = mean(fp_ratio0{ir});
                temp_op_mean(ir) = mean(op_ratio0{ir});
                temp_ofp_mean(ir) = mean(ofp_ratio0{ir});
                temp_p_std(ir) = std0(p_ratio0{ir});
                temp_fp_std(ir) = std0(fp_ratio0{ir});
                temp_op_std(ir) = std0(op_ratio0{ir});
                temp_ofp_std(ir) = std0(ofp_ratio0{ir});
            end

                mean_p   = cat(1,mean_p,temp_p_mean);
                mean_fp  = cat(1,mean_fp,temp_fp_mean);
                mean_op  = cat(1,mean_op,temp_op_mean);
                mean_ofp = cat(1,mean_ofp,temp_ofp_mean);
                std_p   = cat(1,std_p,temp_p_std);
                std_fp  = cat(1,std_fp,temp_fp_std);
                std_op  = cat(1,std_op,temp_op_std);
                std_ofp = cat(1,std_ofp,temp_ofp_std);
            hist_text = cat(2,hist_text,{[image_folder,' foci'],[image_folder,' fake']});
            
            clear foci_data fake_data p_ratio0 fp_ratio0 op_ratio0 ofp_ratio0 h
            
            data73a = cat(1,data73a,out73a);
            data73b = cat(1,data73b,out73b);
            data74  = cat(1,data74 ,out74 );
            data75  = cat(1,data75 ,out75 );
            data76a = cat(1,data76a,out76a);
            data76b = cat(1,data76b,out76b);
            data77  = cat(1,data77 ,out77 );
            data_legend0 = cat(2,data_legend0,legend_text0);
            data_legend  = cat(2,data_legend ,legend_text );
            data_legend1 = cat(2,data_legend1,legend_text1);
        end
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Histogram plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hist_foci_I = hist_I+[1:2:(2*size(mean_p,1)-1)];
    hist_fake_I = hist_I+[2:2:(2*size(mean_p,1))];
    for ir = 1:length(r_size)
        figure(78)
        subplot(2,ceil(length(r_size)/2),ir)
            hold on
            bar(hist_foci_I,mean_p(:,ir)*100,'BarWidth',0.4,'FaceColor',ccode_foci{list_G})
            errorbar(hist_foci_I,mean_p(:,ir)*100,std_p(:,ir)*100,'color','k','LineStyle','none')
            bar(hist_fake_I,mean_fp(:,ir)*100,'BarWidth',0.4,'FaceColor',ccode_fake{list_G})
            errorbar(hist_fake_I,mean_fp(:,ir)*100,std_fp(:,ir)*100,'color','k','LineStyle','none')
        figure(79)
        subplot(2,ceil(length(r_size)/2),ir)
            hold on
            bar(hist_foci_I,mean_op(:,ir),'BarWidth',0.4,'FaceColor',ccode_foci{list_G})
            errorbar(hist_foci_I,mean_op(:,ir),std_op(:,ir),'color','k','LineStyle','none')
            bar(hist_fake_I,mean_ofp(:,ir),'BarWidth',0.4,'FaceColor',ccode_fake{list_G})
            errorbar(hist_fake_I,mean_ofp(:,ir),std_ofp(:,ir),'color','k','LineStyle','none')
    end        
    hist_I = hist_fake_I(end);
        
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Averaging: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    enrichment_mean(data73a,data73b,data74,data75,data76a,data76b,data77,foci_data0,fake_data0,p_ratio00,fp_ratio00,op_ratio00,ofp_ratio00,data_legend0,data_legend,data_legend1,group_name,t_absolute,r_size);    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% Output: %%%============================================================
    result_folder = [output_folder,group_name,'/'];
    if exist(result_folder) ~= 7
        mkdir(result_folder);
    end
    saveas(73,[result_folder,group_name,protein_add,local_add,foci_add,figure_tail]);
    saveas(74,[result_folder,group_name,protein_add,local_add,hist_add,figure_tail]);
    saveas(75,[result_folder,group_name,protein_add,local_add,RNA_add,figure_tail]);
    saveas(76,[result_folder,group_name,protein_add,local_add,abs_add,foci_add,figure_tail]);
    saveas(77,[result_folder,group_name,protein_add,local_add,abs_add,hist_add,figure_tail]);
    close(73)
    close(74)
    close(75)
    close(76)
    close(77)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    hist_legend = cat(2,hist_legend,{[group_name,' foci mean'],[group_name,' foci std'],[group_name,' fake mean'],[group_name,' fake std']});

end



if t_absolute
    unit1 = '#';
    unit2 = 'M';
else
    unit1 = 'A.U.';
    unit2 = 'A.U.';
end

for ir = 1:length(r_size)
    figure(78)
    subplot(2,ceil(length(r_size)/2),ir)
        hold on
        title(['Mean enrichment intensity histogram (r = ',num2str(r_size(ir)),') '])
        ylabel('Mean enrichment (%)')
%         xticklabel_rotate([1:hist_I],90,hist_text,'Interpreter','none');
        legend(hist_legend)
        legend('off')
    saveas(78,[output_folder,protein_add,local_add,hist_add,figure_tail]);
    figure(79)
    subplot(2,ceil(length(r_size)/2),ir)
        hold on
        title(['Mean enrichment # histogram (r = ',num2str(r_size(ir)),') '])
        ylabel(['Mean enrichment (',unit1,')'])
%         xticklabel_rotate([1:hist_I],90,hist_text,'Interpreter','none');
        legend(hist_legend)
        legend('off')
    saveas(79,[output_folder,protein_add,local_add,abs_add,hist_add,figure_tail]);
end        
    










%toc