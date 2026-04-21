close all
clear all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_folder = '\enrichment_result\OreR_FISHIF_hb48TMR_BcdA488_60X_2D_cycle11-12_000-100_bgCon2_nonorm_all_newscaling\r=4\';
data_folder = '\enrichment_result\9491_FISHIF_lacZTMR_nulloA647_BcdA488_60X_2D_cycle12-13_000-100_bgCon2_nonorm_all\r=4\';
f_list = dir(data_folder);
true_list = find([f_list.isdir]);
sub_list = {f_list(true_list(3:end)).name};
data_tail = '*.mat';
output_tail1 = '.fig';
output_tail2 = '.eps';
output_tail3 = '.png';

sub_pos0 = [3,9];
sub_pos = [3,3];
N_bin = 25;

Cw = 1;
cbin = 1;
C_bin = Cw:cbin:40-Cw;
C_min = C_bin-Cw;
C_max = C_bin+Cw;
Cw2 = 5;
cbin2 = 5;
C_bin2 = Cw:cbin2:60-Cw;
C_min2 = C_bin2-Cw2;
C_max2 = C_bin2+Cw2;
C_renorm = 1e-9;

EL_min = 0.25;
EL_max = 0.7;
EL_range = [num2str(EL_min*100,'%03u'),'-',num2str(EL_max*100,'%03u')];
out_name1 = ['enrich_all_EL',EL_range];
out_name2 = ['var_all_EL',EL_range];
out_name4 = ['TX_all_EL',EL_range];
out_name3 = ['mean_all_EL',EL_range];

unit1 = '#';
unit2 = 'nM';
real_name = 'hb';
control_name = 'nullo';
pro_name = 'Bcd';
color_code = mycolors(3);
legend_txt = {'foci-bg','fake-bg','foci-fake'};
I_switch = [0,0,1];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure;   %%% figure initiation for individual mean plots
set(h1,'Units','inches')
set(h1,'Position',[2,0.5,18,10])
h2 = figure;   %%% figure initiation for individual var plots
set(h2,'Units','inches')
set(h2,'Position',[2,0.5,18,10])
h4 = figure;   %%% figure initiation for individual TX plots
set(h4,'Units','inches')
set(h4,'Position',[2,0.5,18,10])
h3 = figure;   %%% figure initiation for average plots
set(h3,'Units','inches')
set(h3,'Position',[4,0.5,8.5,11])
ha1 = zeros(0);
ha2 = zeros(0);
ha4 = zeros(0);
ha31 = zeros(0);
ha32 = zeros(0);
ha33 = zeros(0);

I_sample = 0;
for I_sub = 1:length(sub_list)
    %%% Data loading:
    sub_name = sub_list{I_sub};
    data_file = dir([data_folder,sub_name,'\',data_tail]);
    fname = data_file.name;
    load([data_folder,sub_name,'\',fname],'enrich_foci','var_foci','TX_foci','folder_all')
if ~isempty(enrich_foci{2})
    enrich_raw = enrich_foci{2};
    var_raw = var_foci{2};
    TX_raw = TX_foci{2};
    Itrue = (enrich_raw(:,end) >= EL_min) & (enrich_raw(:,end) <= EL_max) & (enrich_raw(:,1) > 0);
    enrich_raw = enrich_raw(Itrue,1:end-1);
    enrich_raw(:,1:2:end) = enrich_raw(:,1:2:end)/C_renorm;
    var_raw = var_raw(Itrue,1:end-1);
    var_raw(:,1:3:end) = var_raw(:,1:3:end)/C_renorm;
    TX_raw = TX_raw(Itrue,:);
    TX_raw(:,3:3:end) = TX_raw(:,3:3:end)/C_renorm;

    %%% Data segmentation and individual plotting:
    N_seg = [1;cumsum(diff(enrich_raw(:,1)) < 0)+1];
    for J_seg = 1:max(N_seg)
        I_sample = I_sample+1;
        folder_name = [folder_all{J_seg}(1:find(folder_all{J_seg} == '/' | folder_all{J_seg} == '\',1)),folder_all{J_seg}(end-7:end)];
        figure(h1)
        ha1(I_sample) = subaxis(sub_pos0(1),sub_pos0(2),J_seg,I_sub, 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        for I_plot = 1:size(enrich_raw,2)/2-1
            enrich_out = enrich_raw(N_seg == J_seg,(I_plot*2-1):(I_plot*2));
            [xbin,ybin,xerr,yerr] = equal_dist(enrich_out(:,1),enrich_out(:,2),C_min,C_max);
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2,'DisplayName',legend_txt{I_plot})
            hold on
        end
%         title(['Sample ',num2str(I_sample),', Cycle ',sub_name])
        title([folder_name,char(10),'Cycle ',sub_name],'Interpreter','none')
        ylabel([pro_name,' enrichment (',unit1,')'])
        xlabel(['[',pro_name,'] (',unit2,')'])
        legend(legend_txt(1:I_plot))
        legend('off')

        figure(h2)
        ha2(I_sample) = subaxis(sub_pos0(1),sub_pos0(2),J_seg,I_sub, 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        for I_plot = 1:size(var_raw,2)/3-1
            var_out = var_raw(N_seg == J_seg,(I_plot*3-2):(I_plot*3));
            [xbin,ybin,xerr,yerr] = equal_dist_var(var_out(:,1),var_out(:,2),var_out(:,3),I_switch(I_plot),C_min,C_max);
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2,'DisplayName',legend_txt{I_plot})
            hold on
        end
%         title(['Sample ',num2str(I_sample),', Cycle ',sub_name])
        title([folder_name,char(10),'Cycle ',sub_name],'Interpreter','none')
        ylabel([pro_name,' enrichment variance (',unit1,')'])
        xlabel(['[',pro_name,'] (',unit2,')'])
        legend(legend_txt(1:I_plot))
        legend('off')

        figure(h4)
        ha4(I_sample) = subaxis(sub_pos0(1),sub_pos0(2),J_seg,I_sub, 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        for I_plot = 1:size(TX_raw,2)/3
            TX_out = TX_raw(N_seg == J_seg,(I_plot*3-2):(I_plot*3));
            [ybin,xbin,yerr,xerr] = equal_dist(TX_out(:,2),TX_out(:,1),C_min2,C_max2);
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2,'DisplayName',legend_txt{I_plot})
            hold on
        end
%         title(['Sample ',num2str(I_sample),', Cycle ',sub_name])
        title([folder_name,char(10),'Cycle ',sub_name],'Interpreter','none')
        xlabel([pro_name,' enrichment variance (',unit1,')'])
        ylabel([real_name,' TX level (',unit1,')'])
        legend(legend_txt(1:I_plot))
        legend('off')
    end
    
    %%% Average data plotting:
    figure(h3)
    ha31(I_sub) = subaxis(sub_pos(1),sub_pos(2),1,I_sub, 'Spacing', 0, 'PaddingRight',0.04, 'PaddingLeft',0.04, 'PaddingTop',0.035, 'PaddingBottom',0.035, 'MarginLeft', 0.04, 'MarginRight', 0.005, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        for I_plot = 1:size(enrich_raw,2)/2-1
            enrich_out = enrich_raw(:,(I_plot*2-1):(I_plot*2));
            [xbin,ybin,xerr,yerr] = equal_dist(enrich_out(:,1),enrich_out(:,2),C_min,C_max);
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2,'DisplayName',legend_txt{I_plot})
            hold on
        end
        title(['Cycle ',sub_name],'FontName','Arial','FontSize',10,'FontWeight','bold')
        ylabel([pro_name,' enrichment (',unit1,')'],'FontName','Arial','FontSize',10,'FontWeight','bold')
        xlabel(['[',pro_name,'] (',unit2,')'],'FontName','Arial','FontSize',10,'FontWeight','bold')
        set(gca,'FontName','Arial','FontSize',10)
        
    ha32(I_sub) = subaxis(sub_pos(1),sub_pos(2),2,I_sub, 'Spacing', 0, 'PaddingRight',0.04, 'PaddingLeft',0.04, 'PaddingTop',0.035, 'PaddingBottom',0.035, 'MarginLeft', 0.04, 'MarginRight', 0.005, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        for I_plot = 1:size(var_raw,2)/3-1
            var_out = var_raw(:,(I_plot*3-2):(I_plot*3));
            [xbin,ybin,xerr,yerr] = equal_dist_var(var_out(:,1),var_out(:,2),var_out(:,3),I_switch(I_plot),C_min,C_max);
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2,'DisplayName',legend_txt{I_plot})
            hold on
        end
        title(['Cycle ',sub_name],'FontName','Arial','FontSize',10,'FontWeight','bold')
        ylabel([pro_name,' enrichment variance (',unit1,')'],'FontName','Arial','FontSize',10,'FontWeight','bold')
        xlabel(['[',pro_name,'] (',unit2,')'],'FontName','Arial','FontSize',10,'FontWeight','bold')
        set(gca,'FontName','Arial','FontSize',10)
        
    ha33(I_sub) = subaxis(sub_pos(1),sub_pos(2),3,I_sub, 'Spacing', 0, 'PaddingRight',0.04, 'PaddingLeft',0.04, 'PaddingTop',0.035, 'PaddingBottom',0.035, 'MarginLeft', 0.04, 'MarginRight', 0.005, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        for I_plot = 1:size(TX_raw,2)/3
            TX_out = TX_raw(:,(I_plot*3-2):(I_plot*3));
            [ybin,xbin,yerr,xerr] = equal_dist(TX_out(:,2),TX_out(:,1),C_min2,C_max2);
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,color_code(I_plot,:),color_code(I_plot,:),'o','LineWidth',2,'DisplayName',legend_txt{I_plot})
            hold on
        end
        title(['Cycle ',sub_name],'FontName','Arial','FontSize',10,'FontWeight','bold')
        xlabel([pro_name,' enrichment variance (',unit1,')'],'FontName','Arial','FontSize',10,'FontWeight','bold')
        ylabel([real_name,' TX level (',unit1,')'],'FontName','Arial','FontSize',10,'FontWeight','bold')
        set(gca,'FontName','Arial','FontSize',10)
end
end
            
linkaxes(ha1)
    axes(ha1(1))
    ylim([-0.2,0.4])
linkaxes(ha2)
    axes(ha2(1))
    ylim([-0.2,0.2])
linkaxes(ha4)
%     axes(ha4(1))
%     ylim([-0.2,0.2])
linkaxes(ha31)
    axes(ha31(1))
    ylim([-0.2,0.4])
linkaxes(ha32)
    axes(ha32(1))
    ylim([-0.2,0.2])
linkaxes(ha33)
%     axes(ha33(1))
%     ylim([-0.2,0.2])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(h1,[data_folder,out_name1,output_tail1])
set(h1,'PaperPositionMode','auto')
export_fig([data_folder,out_name1,output_tail2],'-nocrop','-transparent',h1)
export_fig([data_folder,out_name1,output_tail3],'-nocrop','-transparent',h1)
            
saveas(h2,[data_folder,out_name2,output_tail1])
set(h2,'PaperPositionMode','auto')
export_fig([data_folder,out_name2,output_tail2],'-nocrop','-transparent',h2)
export_fig([data_folder,out_name2,output_tail3],'-nocrop','-transparent',h2)
            
saveas(h4,[data_folder,out_name4,output_tail1])
set(h4,'PaperPositionMode','auto')
export_fig([data_folder,out_name4,output_tail2],'-nocrop','-transparent',h4)
export_fig([data_folder,out_name4,output_tail3],'-nocrop','-transparent',h4)
            
saveas(h3,[data_folder,out_name3,output_tail1])
set(h3,'PaperPositionMode','auto')
export_fig([data_folder,out_name3,output_tail2],'-nocrop','-transparent',h3)
export_fig([data_folder,out_name3,output_tail3],'-nocrop','-transparent',h3)
            
            
            
            
            
            
            
            