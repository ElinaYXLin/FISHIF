function TXnoise_plot(nucleus_RNA_profile,N_cycle,image_folder,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TXnoise_plot: A function to plot the noise level of transcription profile
%% 1. Embryo binning based on EL;
%% 2. Calculate and plot sigma/mean vs 1/mean;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nucleus_RNA_profile (input): nucleus RNA profile(:,[1,4]) = [EL,TX];
%% N_cycle (input): cycle of embryo;
%% image_folder (input): image folder;
%% varargin{1} (input): output figure number;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin{1})
    h1 = 101;
    h2 = 102;
else
    h1 = varargin{1}(1);
    h2 = varargin{1}(2);
end

EL_min0 = 0.1;
EL_max0 = 0.5;
EL_bin = 0.05;
EL_delta = 0.025;
EL_min = EL_min0:EL_delta:(EL_max0-EL_bin);
EL_max = (EL_min0+EL_bin):EL_delta:EL_max0;
meanTX = zeros(size(EL_min));
stdTX = zeros(size(EL_min));

ccode = 'crgbky';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate the noise level: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:length(EL_min)
    Itrue = (nucleus_RNA_profile(:,1) >= EL_min(ii)) & (nucleus_RNA_profile(:,1) <= EL_max(ii));
    meanTX(ii) = mean(nucleus_RNA_profile(Itrue,4));
    stdTX(ii) = std(nucleus_RNA_profile(Itrue,4));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
icolor = N_cycle-9;
icolor(icolor < 1) = 1;icolor(icolor > 6) = 6;
figure(h1)
clf
subplot(1,2,1)
    plot(meanTX,stdTX./meanTX,'Color',ccode(icolor),'Marker','.','MarkerSize',15,'LineStyle','none','DisplayName',image_folder)
    set(gca,'FontName','Arial','FontSize',16)
%     ylim([0,1])
    xlabel('mean','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel('Noise','FontName','Arial','FontSize',16,'FontWeight','bold')
    title([image_folder,': Cycle ',num2str(N_cycle)],'FontName','Arial','FontSize',16,'FontWeight','bold','Interpreter','none')
subplot(1,2,2)
    plot(1./meanTX,(stdTX./meanTX).^2,'Color',ccode(icolor),'Marker','.','MarkerSize',15,'LineStyle','none','DisplayName',image_folder)
    set(gca,'FontName','Arial','FontSize',16)
%     ylim([0,1])
    xlabel('mean^-^1','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel('Noise^2','FontName','Arial','FontSize',16,'FontWeight','bold')
%     title('Noise level vs mean^-^1 (cyan: C10, red: C11, green: C12, blue: C13, black: C14)','FontName','Arial','FontSize',16,'FontWeight','bold')

figure(h2)
subplot(1,2,1)
    plot(meanTX,stdTX./meanTX,'Color',ccode(icolor),'Marker','.','MarkerSize',15,'LineStyle','none','DisplayName',image_folder)
    set(gca,'FontName','Arial','FontSize',16)
    hold on
    ylim([0,1])
    xlabel('mean','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel('Noise','FontName','Arial','FontSize',16,'FontWeight','bold')
    title('Noise level vs mean (cyan: C10, red: C11, green: C12, blue: C13, black: C14)','FontName','Arial','FontSize',16,'FontWeight','bold')
subplot(1,2,2)
    plot(1./meanTX,(stdTX./meanTX).^2,'Color',ccode(icolor),'Marker','.','MarkerSize',15,'LineStyle','none','DisplayName',image_folder)
    set(gca,'FontName','Arial','FontSize',16)
    hold on
    ylim([0,1])
    xlabel('mean^-^1','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel('Noise^2','FontName','Arial','FontSize',16,'FontWeight','bold')
%     title('Noise level vs mean^-^1 (cyan: C10, red: C11, green: C12, blue: C13, black: C14)','FontName','Arial','FontSize',16,'FontWeight','bold')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

