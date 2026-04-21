clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to print embryo analysis results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
% out_folder = 'Results_protein2/';
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
% signal2_add = '_nullo';
% signal2_add = '_hb-y';
protein_add = '_protein';
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
reg_add = '_regulation';
enrich_add = '_enrich';
binN_add = '_binN';
binC_add = '_binC';
para_add = '_para';
prefit_add = '_prefit';

embryo_type = {'T'};
embryo_name = {'Cad-Bcd-hb-gal4-16249-1-5M'};
sub_ij = [6,7];
cycle_range = 10:13;
rsize0 = 2;
nM_con = 1;
Nr = 20;

bin_min = 0:0.025:0.95;   %%% bin_min for EL
bin_max = 0.05:0.025:1;   %%% bin_max for EL
bin_min2 = 0:2.5:100;   %%% bin_min for [Bcd]
bin_max2 = 5:2.5:105;   %%% bin_max for [Bcd]
% % bin_min2 = 0:0.2:10;   %%% bin_min for [Bcd]
% % bin_max2 = 0.4:0.2:10.4;   %%% bin_max for [Bcd]
EL_range = [0.25,0.75];
EL_range2 = [0,1];
z_range = [2,22];
Nbin = 50;
rover = 0.6;
% NbinE = 20;
roverE = 0.6;
% NbinE = 20;
% roverE = 0.6;
NbinE = 30;   %%% Equal population binning for enrichment
NmoveE = 250;
RNAmax = 150;
xlim0 = [0,1]; xlim_bin = 0.01;
ylim2 = [0,100];
ylim2b = [0,200]; 
% ylim3 = [0,100];
ylim3 = [0,30];

fsize = 8;

posf = [1,0.5,8.5,11];

pos7 = [0.75,0.60,4.5,0.75];
pos6 = [0.75,1.35,4.5,1.90]; pos6b = [5.35,1.35,2.75,1.90];
pos5 = [0.75,3.30,4.5,0.75];
pos4 = [0.75,4.05,4.5,1.90]; pos4b = [5.35,4.05,2.75,1.90];
pos3 = [0.75,6.00,4.5,0.75];
pos2 = [0.75,6.75,4.5,1.90]; pos2b = [5.35,6.75,0.25,1.90];
pos1 = [0.75,8.65,4.5,1.90];


% % posf = [1,0.5,11,8.5];
% % 
% % pos7 = [0.5,0.60,4.5,0.75];
% % pos5 = [5.5,0.60,4.5,0.75];
% % pos6 = [0.5,1.35,4.5,2.00]; pos6b = [5.05,1.35,2.75,1.90];
% % pos4 = [5.5,1.35,4.5,2.00]; pos4b = [5.6,4.05,2.75,1.90];
% % pos3 = [0.5,3.45,4.5,0.75];
% % pos32 = [5.5,3.45,4.5,0.75];
% % pos2 = [0.5,4.20,4.5,2.00]; pos2b = [5.6,6.75,0.25,1.90];
% % pos22 = [5.5,4.20,4.5,2.00]; pos22b = [5.6,6.75,0.25,1.90];
% % pos1 = [5.5,6.20,4.5,2.00];

flip_axis = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for  iem = 1:length(embryo_type)
    
%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi(embryo_type{iem},folder_list(:,6)),:);
    [N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for list_I = 1:N1
        [sub_num, sub_list, sub_all] = xlsread([folder_list{list_I,1},in_folder,input_name]);
        channel_name = eval(folder_list{list_I,5});
        if size(sub_all,2) >= 17
            I_se = ~strcmp('F',sub_all(:,17));
    %         I_se = strcmp('T',sub_all(:,17));
            sub_num = sub_num(I_se,:);
            sub_list = sub_list(I_se,:);
        end
        [M1,M2] = size(sub_list);

        for list_J = 1:M1
            image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
            result_folder = [folder_list{list_I,1},out_folder];
            load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_protein_profile','nucleus_protein_profile_ab','nucleus_RNA_profile','foci_RNA_profile','nucleus_signal2_profile','foci_signal2_profile','max_image','DAPI_channel','protein_channel');
            nucleus_protein_profile_ab(:,2) = nucleus_protein_profile_ab(:,2)/1e-9;
            foci_RNA_profile = foci_RNA_profile(foci_RNA_profile(:,2) > 0,:);
            Iz = nucleus_protein_profile(:,4) >= z_range(1) & nucleus_protein_profile(:,4) <= z_range(2);
            N_cycle = sub_num(list_J,13);   %%% Calculate the nuclear cycle number
            title0 = [folder_list{list_I,1}(find(folder_list{list_I,1}(1:end-1) == '\' | folder_list{list_I,1}(1:end-1) == '/',1,'last')+1:end),sub_list{list_J,3}];
            
            
%% Figure initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure(1); set(1,'Units','Inches','Position',posf,'PaperPositionMode', 'auto')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% DAPI image replot: %%%=================================================
            figure(1)
            ha10 = axes;
            set(ha10,'Units','Inches','Position',pos1)
%                 im0 = double(max_image(:,:,DAPI_channel));
%                 im0 = max(double(max_image(:,:,[2,5])),[],3);
                im0a = double(max_image(:,:,2));
                im0b = double(max_image(:,:,5));
                im0 = (im0a/max(im0a(:))+im0b/max(im0b(:)))/2;
                imshow(repmat(im0,[1,1,3])*4/max(im0(:)))
                title([title0,': Cycle ',num2str(N_cycle)],'FontSize',fsize,'Interpreter','none')
                cmap0 = colormap;
            set(ha10,'Units','Inches','Position',pos1)
            
%%% Bcd concentration analysis: %%%========================================
            %%% Bcd fate plot
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos2)
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
                cmap = colormap;
                ha0 = get(h0,'Children');
                copyaxes(ha0(end),ha1,true);
            figure(1)
            ha2 = axes;
            set(ha2,'Units','Inches','Position',pos2b);
                copyaxes(ha0(1),ha2,true);
                colormap(cmap);
%                 colormap(ha10,cmap0);
                ylabel(out_folder)
                close(h0)
            set(ha1,'Units','Inches','Position',pos2)

            %%% [Bcd] vs. EL
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos3)
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_protein_profile_ab(Iz,1),nucleus_protein_profile_ab(Iz,2),bin_min,bin_max);
                Itrue = (xbin >= EL_range2(1)) & (xbin <= EL_range2(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                plot(nucleus_protein_profile_ab(Iz,1),nucleus_protein_profile_ab(Iz,2),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                xlim(xlim0)
                ylim(ylim3)
%                 xlabel('EL')
                ylabel('[Protein] (nM)','FontSize',fsize)
                hold off
            set(ha1,'Units','Inches','Position',pos3,'XTick',[],'FontSize',fsize)

%%% Foci number analysis: %%%==============================================
            %%% Dostatni plot for RNA channel
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos4)
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0,ha1,true);
                close(h0)
            set(ha1,'Units','Inches','Position',pos4)

            %%% Dostatni plot for RNA2 channel
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos6)
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0,ha1,true);
                close(h0)
            set(ha1,'Units','Inches','Position',pos6)

% % % %             %%% #foci/nu distribution for RNA channel
% % % %             figure(1)
% % % %             ha1 = axes;
% % % %             set(ha1,'Units','Inches','Position',pos5)
% % % %                 h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
% % % %                 ha0 = get(h0,'Children');
% % % %                 copyaxes(ha0(2),ha1,true);
% % % %                 close(h0)
% % % %             set(ha1,'Units','Inches','Position',pos5,'XTick',[])
% % % % 
% % % %             %%% #foci/nu distribution for RNA2 channel
% % % %             figure(1)
% % % %             ha1 = axes;
% % % %             set(ha1,'Units','Inches','Position',pos7)
% % % %                 h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,num_add,D3_add,figure_tail]);
% % % %                 ha0 = get(h0,'Children');
% % % %                 copyaxes(ha0(2),ha1,true);
% % % %                 close(h0)
% % % %             set(ha1,'Units','Inches','Position',pos7,'XTick',[])

% %             %%% #foci/nu vs. EL for RNA channel
% %             figure(1)
% %             ha1 = axes;
% %             set(ha1,'Units','Inches','Position',pos5)
% %                 [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_RNA_profile(Iz,1),nucleus_RNA_profile(Iz,3),bin_min,bin_max);
% %                 Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
% %                 xbin0 = xbin(Itrue);
% %                 ybin0 = ybin(Itrue);
% %                 if isempty(foci_prefit)
% %                     [rmax,I_max] = max(ybin0);
% %                     beta0 = [hc0,rmax,0];
% %                     try
% %                         [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
% %                     catch err
% %                         beta1 = nan(size(beta0));
% %                     end
% %                 else
% %                     beta1 = foci_prefit(i00,2:end);
% %                 end
% %                 beta_all1 = [beta_all1;beta1];
% %                 errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
% %                 hold on
% %                 plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta1,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
% %                 text(mean(xlim0),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
% %                 xlim(xlim0)
% %                 ylim(ylim1)
% %                 xlabel('EL')
% %                 ylabel('#foci/nu')
% %                 hold off
% %             set(ha1,'Units','Inches','Position',pos5,'XTick',[])
% % 
% %             %%% #foci/nu vs. EL for RNA2 channel
% %             figure(1)
% %             ha1 = axes;
% %             set(ha1,'Units','Inches','Position',pos7)
% %                 [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_signal2_profile(Iz,1),nucleus_signal2_profile(Iz,3),bin_min,bin_max);
% %                 Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
% %                 xbin0 = xbin(Itrue);
% %                 ybin0 = ybin(Itrue);
% %                 if isempty(foci2_prefit)
% %                     [rmax,I_max] = max(ybin0);
% %                     beta0 = [hc0,rmax,0];
% %                     try
% %                         [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
% %                     catch err
% %                         beta1 = nan(size(beta0));
% %                     end
% %                 else
% %                     beta1 = foci2_prefit(i00,2:end);
% %                 end
% %                 beta_all1b = [beta_all1b;beta1];
% %                 errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
% %                 hold on
% %                 plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta1,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
% %                 text(mean(xlim0),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
% %                 xlim(xlim0)
% %                 ylim(ylim1)
% %                 xlabel('EL')
% %                 ylabel('#foci/nu')
% %             set(ha1,'Units','Inches','Position',pos7,'XTick',[])

            %%% #RNA/nu vs. EL for RNA channel
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos5)
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_RNA_profile(Iz,1),nucleus_RNA_profile(Iz,4),bin_min,bin_max);
                Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                plot(nucleus_RNA_profile(Iz,1),nucleus_RNA_profile(Iz,4),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                xlim(xlim0)
                ylim(ylim2)
%                 xlabel('EL')
                ylabel('#RNA/nu','FontSize',fsize)
                hold off
            set(ha1,'Units','Inches','Position',pos5,'XTick',[],'FontSize',fsize)

            %%% #RNA/nu vs. EL for RNA2 channel
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos7)
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_signal2_profile(Iz,1),nucleus_signal2_profile(Iz,4),bin_min,bin_max);
                Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                plot(nucleus_signal2_profile(Iz,1),nucleus_signal2_profile(Iz,4),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                xlim(xlim0)
                ylim(ylim2b)
                xlabel('EL','FontSize',fsize)
                ylabel('#RNA/nu','FontSize',fsize)
                hold off
            set(ha1,'Units','Inches','Position',pos7,'FontSize',fsize)

            

%%% Foci number analysis: %%%==============================================
            xmin = 0.30;
            xmax = 0.40;
            ymin = 0.45;
            ymax = 0.55;
            imsize = size(max_image);
            
            xmin1 = 0.25;
            xmax1 = 0.35;
            
            flip_EL = max_image(round(imsize(1)*ymin):round(imsize(1)*ymax),round(imsize(2)*xmin):round(imsize(2)*xmax),protein_channel) < max_image(round(imsize(1)*ymin):round(imsize(1)*ymax),round(imsize(2)*(1-xmax)):round(imsize(2)*(1-xmin)),protein_channel);
            flip_all = xor(flip_axis,flip_EL);
            
            if flip_all
                xlim00 = [1-xmax1,1-xmin1];
            else
                xlim00 = [xmin1,xmax1];
            end
            xlimA = round(xlim00*imsize(2));
            ylimA = round([imsize(1)/2-(xlimA(2)-xlimA(1))*pos4b(4)/pos4b(3)/2,imsize(1)/2+(xlimA(2)-xlimA(1))*pos4b(4)/pos4b(3)/2]);
                

            %%% RNA1 foci recognition plots
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos4b)
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                data0 = get(get(ha0(1),'Children'),'CData');
            axes(ha1)
                imshow(data0(ylimA(1):ylimA(2),xlimA(1):xlimA(2),:)*4)
%             xlim(xlimA)
%             ylim(ylimA)
%             axis equal
% %             xlabel(RNA_add)
            close (h0)

            %%% RNA2 foci recognition plots
            figure(1)
            ha1 = axes;
            set(ha1,'Units','Inches','Position',pos6b)
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,seg_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                data0 = get(get(ha0(1),'Children'),'CData');
            axes(ha1)
                imshow(data0(ylimA(1):ylimA(2),xlimA(1):xlimA(2),:)*8)
%             xlim(xlimA)
%             ylim(ylimA)
%             axis equal
% %             xlabel(RNA_add)
            close (h0)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

%% Print figure: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(1,'Renderer','OpenGL')
% %             orient landscape
            print
            close(1)
        end
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end




