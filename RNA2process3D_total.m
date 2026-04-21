clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to summarize the result of dual-FISH labeling experiment %%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNA2list.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'Results/';

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
all_add = '_all';

embryo_type = {'m1','m2','m3','m4'};
embryo_name = {'16249-1-1M','16249-1-5M','16249-2-1M','16249-2-2M'};
sub_ij = [3,4];

bin_min = 0:0.025:0.95;
bin_max = 0.05:0.025:1;
reg_range = [0.25,0.75];
xlim0 = [0,1]; xlim_bin = 0.01;
ylim1 = [0,3];
ylim2 = [0,200];
ylim3 = [0,100];
hc0 = [-6,0.5];

hist_binxy = 0:0.1:10;
hist_binz = -5:5;

ana_folder = 'RNA2process3D_result/';
xls_title = {'Name','Cycle','h','EL0','ymax','ymin'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for  iem = 1:length(embryo_type)
    
%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi(embryo_type{iem},folder_list(:,6)),:);
    [N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    figure(1); maximize(1); set(1,'name',['Dostatni plot (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(2); maximize(2); set(2,'name',['Dostatni plot (Marker): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(3); maximize(3); set(3,'name',['#foci/nu distribution (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(4); maximize(4); set(4,'name',['#foci/nu distribution (Marker): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(5); maximize(5); set(5,'name',['#foci/nu (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(6); maximize(6); set(6,'name',['#foci/nu (Marker): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(7); maximize(7); set(7,'name',['#RNA/nu: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(8); maximize(8); set(8,'name',['#RNA/marker: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(9); maximize(9); set(9,'name',['dist_xy: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(10); maximize(10); set(10,'name',['dist_z: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(11); set(11,'name',['foci-marker match: ',embryo_name{iem}],'Units','inches','Position',[1,1,10,5]); set(11,'Units','normal','PaperPositionMode', 'auto')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i00 = 0;
    all_foci_xy = zeros(0);
    all_foci_z = zeros(0);
    beta_all1 = zeros(0);
    beta_all2 = zeros(0);
    beta_all3 = zeros(0);
    name_all = cell(0);
    cycle_all = zeros(0);

    for list_I = 1:N1
        [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
        channel_name = eval(folder_list{list_I,5});
        [M1,M2] = size(sub_list);

        for list_J = 1:M1
            i00 = i00+1;
            image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
            result_folder = [folder_list{list_I,1},out_folder];
            load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_RNA_profile','foci_RNA_profile','nucleus_signal2_profile','foci_signal2_profile','foci_mismatch','dxy0','marker_RNA_profile');
            N_cycle = sub_num(list_J,13);   %%% Calculate the nuclear cycle number
            all_foci_xy = [all_foci_xy;dxy0];
            all_foci_z = [all_foci_z;foci_mismatch(:,7)];
            name_all = cat(1,name_all,{image_folder});
            cycle_all = [cycle_all;N_cycle];

%%% Foci number analysis: %%%==============================================
            %%% Dostatni plot for RNA channel
            figure(1)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0,ha1,true);
                close(h0)
                axes(ha1);
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')

            %%% Dostatni plot for marker channel
            figure(2)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0,ha1,true);
                close(h0)
                axes(ha1);
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')

            %%% #foci/nu distribution for RNA channel
            figure(3)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0(2),ha1,true);
                close(h0)
                axes(ha1);
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')

            %%% #foci/nu distribution for marker channel
            figure(4)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,num_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0(2),ha1,true);
                close(h0)
                axes(ha1);
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')

            %%% #foci/nu for RNA channel
            figure(5)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,3),bin_min,bin_max);
                Itrue = (xbin >= reg_range(1)) & (xbin <= reg_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                [rmax,I_max] = max(ybin0);
                beta0 = [hc0,rmax,0];
                try
                    [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
                catch err
                    beta1 = nan(size(beta0));
                end
                beta_all1 = [beta_all1;beta1];
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                hold on
                plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta1,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim0),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim1)
                xlabel('EL')
                ylabel('#foci/nu')
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                hold off

            %%% #foci/nu for marker channel
            figure(6)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_signal2_profile(:,1),nucleus_signal2_profile(:,3),bin_min,bin_max);
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                xlim(xlim0)
                ylim(ylim1)
                xlabel('EL')
                ylabel('#foci/nu')
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')

            %%% #RNA/nu
            figure(7)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,4),bin_min,bin_max);
                Itrue = (xbin >= reg_range(1)) & (xbin <= reg_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                [rmax,I_max] = max(ybin0);
                beta0 = [hc0,rmax,0];
                try
                    [beta2,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
                catch err
                    beta2 = nan(size(beta0));
                end
                beta_all2 = [beta_all2;beta2];
                plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,4),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta2,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim0),mean(ylim2),['h = ',num2str(beta2(1),'%5.2f'),char(10),'EL0 = ',num2str(beta2(2),'%5.2f'),char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim2)
                xlabel('EL')
                ylabel('#RNA/nu')
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                hold off

            %%% #RNA/marker
            figure(8)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(marker_RNA_profile(:,1),marker_RNA_profile(:,2),bin_min,bin_max);
                Itrue = (xbin >= reg_range(1)) & (xbin <= reg_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                [rmax,I_max] = max(ybin0);
                beta0 = [hc0,rmax,0];
                try
                    [beta3,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
                catch err
                    beta3 = nan(size(beta0));
                end
                beta_all3 = [beta_all3;beta3];
                plot(marker_RNA_profile(:,1),marker_RNA_profile(:,2),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta3,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim0),mean(ylim3),['h = ',num2str(beta3(1),'%5.2f'),char(10),'EL0 = ',num2str(beta3(2),'%5.2f'),char(10),'ymax = ',num2str(beta3(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim3)
                xlabel('EL')
                ylabel('#RNA/marker')
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                hold off

            %%% Distribution of xy mismatch
            figure(9)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [nxy,dxy] = hist(dxy0,[hist_binxy,2*hist_binxy(end)-hist_binxy(end-1)]);
                ho1 = bar(dxy,nxy,'hist');
                set(ho1,'FaceColor',0.8*ones(1,3))
                xlabel('dxy')
                ylabel('#')
                xlim(hist_binxy([1,end]))
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')

            %%% Distribution of z mismatch
            figure(10)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [nz,dz] = hist(foci_mismatch(:,7),[2*hist_binz(1)-hist_binz(2),hist_binz,2*hist_binz(end)-hist_binz(end-1)]);
                ho1 = bar(dz,nz,'hist');
                set(ho1,'FaceColor',0.8*ones(1,3))
                xlabel('dz')
                ylabel('#')
                xlim(hist_binz([1,end]))
                title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
        end
    end

%%% Pooled distribution of mismatch
    figure(11)
    subplot(1,2,1);
        [nxy,dxy] = hist(all_foci_xy,[hist_binxy,2*hist_binxy(end)-hist_binxy(end-1)]);
        ho1 = bar(dxy,nxy,'hist');
        set(ho1,'FaceColor',0.8*ones(1,3))
        xlabel('dxy')
        ylabel('#')
        xlim(hist_binxy([1,end]))
        title(['xy distance between foci and marker: ',embryo_name{iem}])
    subplot(1,2,2);
        [nz,dz] = hist(all_foci_z,[2*hist_binz(1)-hist_binz(2),hist_binz,2*hist_binz(end)-hist_binz(end-1)]);
        ho1 = bar(dz,nz,'hist');
        set(ho1,'FaceColor',0.8*ones(1,3))
        xlabel('dz')
        ylabel('#')
        xlim(hist_binz([1,end]))
        title(['z distance between foci and marker: ',embryo_name{iem}])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% Data output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist([ana_folder,embryo_name{iem},'/']) ~= 7
        mkdir([ana_folder,embryo_name{iem},'/']);
    end

    saveas(1,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fate_add,D3_add,figure_tail]);
    saveas(2,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,signal2_add,fate_add,D3_add,figure_tail]);
    saveas(3,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,num_add,D3_add,figure_tail]);
    saveas(4,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,signal2_add,num_add,D3_add,figure_tail]);
    saveas(5,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,figure_tail]);
    saveas(6,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,figure_tail]);
    saveas(7,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,figure_tail]);
    saveas(8,[ana_folder,embryo_name{iem},'/',embryo_name{iem},marker_add,fish_add,int_add,D3_add,figure_tail]);
    saveas(9,[ana_folder,embryo_name{iem},'/',embryo_name{iem},mismatch_add,xy_add,figure_tail]);
    saveas(10,[ana_folder,embryo_name{iem},'/',embryo_name{iem},mismatch_add,z_add,figure_tail]);
    saveas(11,[ana_folder,embryo_name{iem},'/',embryo_name{iem},mismatch_add,all_add,figure_tail]);

    saveas(1,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fate_add,D3_add,figure_tail2]);
    saveas(2,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,signal2_add,fate_add,D3_add,figure_tail2]);
    saveas(3,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,num_add,D3_add,figure_tail2]);
    saveas(4,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,signal2_add,num_add,D3_add,figure_tail2]);
    saveas(5,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,figure_tail2]);
    saveas(6,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,figure_tail2]);
    saveas(7,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,figure_tail2]);
    saveas(8,[ana_folder,embryo_name{iem},'/',embryo_name{iem},marker_add,fish_add,int_add,D3_add,figure_tail2]);
    saveas(9,[ana_folder,embryo_name{iem},'/',embryo_name{iem},mismatch_add,xy_add,figure_tail2]);
    saveas(10,[ana_folder,embryo_name{iem},'/',embryo_name{iem},mismatch_add,z_add,figure_tail2]);
    saveas(11,[ana_folder,embryo_name{iem},'/',embryo_name{iem},mismatch_add,all_add,figure_tail2]);

    %%% Save the Hill fitting results to Excel files:
    out_ex0 = cat(1,xls_title,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all1)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all2)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all3)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},marker_add,fish_add,int_add,D3_add,output_tail],out_ex0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

end




