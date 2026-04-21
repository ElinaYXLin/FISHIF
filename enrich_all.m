clear all
close all

total_folder = 'enrichment_result/';
protein_add = '_protein';
local_add = '_local';
abs_add = '_abs';
rad_add = '_rad';
D3_add = '_3D';
mean_add = '_mean';
spot_add = '_spot';
figure_tail = '.fig';

r_range = 0:9;
N_bin = 25;
N_bin2 = 20;
N_bin3 = 10;

% % % % folder_name = 'OreR_FISHIF_hb48TMR_Act5CA647_BcdA488_60X_2D_cycle12-13_00-100_newEL/';
% % % % folder_name = 'OreR_FISHIF_hb48TMR_BcdA488_nulloA647_60X_2D_cycle12-13_00-100_newEL/';
% % % % folder_name = 'OreR_FISHIF_hb48TMR_BcdA488_60X_2D_cycle12-13_00-100_newEL/';
% % % % folder_name = 'Bcd3_FISHIF_lacZTMR_nulloA647_BcdA488_60X_2D_cycle10-13_000-100_newEL/';
% % % % folder_name = 'OreR_FISHIF_hb48TMR_HbA488_60X_2D_cycle12-13_00-100_newEL/';
% % % % folder_name = '10652-1-3M_FISHIF_lacZ72A647_gal4-48TMR_BcdA488_60X_00-100_newEL/';
% % % % folder_name = '10652-1-3M_FISHIF_all_60X_00-100_newEL/';
folder_name = 'BAC9_FISHIF_all_60X_00-100_newEL/';
% % % % folder_name = 'Bcd3_FISHIF_lacZTMR_nulloA647_HbA488_60X_2D_cycle10-13_000-100_newEL/';
enrich_ex0 = cell(size(r_range));
enrich_ex = cell(size(r_range));

for I0 = 1:length(r_range)
    [enrich_ex0{I0},enrich_ex{I0},legend_text] = Dualprocess5_total([folder_name,'r=',num2str(r_range(I0)),'/'],r_range(I0),N_bin);
%     [enrich_ex0{I0},enrich_ex{I0},legend_text] = Dual2process5_total([folder_name,'r=',num2str(r_range(I0)),'/'],r_range(I0),N_bin,N_bin2,N_bin3);
end


%% Plot r dependence:
close all
unit1 = '#';
unit2 = 'Pixel';

d_data = size(enrich_ex0{1},2)/2;
d_data1 = size(enrich_ex{1}{1},2)/2;
d_data2 = size(enrich_ex{1}{2},2)/2;
d_data3 = size(enrich_ex{1}{3},2)/2;
color_code = mycolors(d_data1);
color_code = [color_code;color_code(1:end-1,:)];

mean_en0 = nan(length(r_range),d_data);
std0_en0 = nan(length(r_range),d_data);
mean_en = nan(length(r_range),d_data);
std0_en = nan(length(r_range),d_data);

figure(1);clf
figure(2);clf
maximize(1)
maximize(2)
for I0 = 1:length(r_range)
    mean_en0(I0,:) = mean(enrich_ex0{I0}(:,2:2:2*d_data));
    std0_en0(I0,:) = std0(enrich_ex0{I0}(:,2:2:2*d_data));
    if ~isempty(enrich_ex{I0}{1})
        mean_en(I0,:) = [mean(enrich_ex{I0}{1}(:,2:2:2*d_data1)),mean(enrich_ex{I0}{2}(:,2:2:2*d_data2)),mean(enrich_ex{I0}{3}(:,2:2:2*d_data3))];
        std0_en(I0,:) = [std0(enrich_ex{I0}{1}(:,2:2:2*d_data1)),std0(enrich_ex{I0}{2}(:,2:2:2*d_data2)),std0(enrich_ex{I0}{3}(:,2:2:2*d_data3))];
    end
end

for I_plot = 1:d_data
    figure(1)
        errorbar(r_range,mean_en0(:,I_plot),std0_en0(:,I_plot),'Color',color_code(I_plot,:),'LineWidth',2)
        hold on
    figure(2)
        errorbar(r_range,mean_en(:,I_plot),std0_en(:,I_plot),'Color',color_code(I_plot,:),'LineWidth',2)
        hold on
end

figure(1)
    title([folder_name,': Enrichment vs integration radius (average over ',num2str(size(enrich_ex0{end},1)),' embryos)'],'Interpreter','none','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel(['Mean enrichment (',unit1,')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
    xlabel(['Radius (',unit2,')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
    legend(legend_text)
    set(gca,'FontName','Arial','FontSize',16,'Linewidth',2)
saveas(1,[total_folder,folder_name,folder_name(1:end-1),protein_add,local_add,abs_add,rad_add,D3_add,mean_add,figure_tail]);

figure(2)
    title([folder_name,': Enrichment vs integration radius (average over ',num2str(size(enrich_ex{end}{1},1)),' foci)'],'Interpreter','none','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel(['Mean enrichment (',unit1,')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
    xlabel(['Radius (',unit2,')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
    legend(legend_text)
    set(gca,'FontName','Arial','FontSize',16,'Linewidth',2)
saveas(2,[total_folder,folder_name,folder_name(1:end-1),protein_add,local_add,abs_add,rad_add,D3_add,mean_add,spot_add,figure_tail]);


