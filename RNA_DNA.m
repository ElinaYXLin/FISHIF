function RNA_DNA(DNA_profile,RNA_profile,image_folder,N_cycle,Nfig)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to plot FISH signal against DNA signal %%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DNA_profile: DNA intensity profile;
%% RNA_profile: RNA foci number, RNA intensity profile;
%% image_folder: input image folder;
%% N_cycle: embryo cycle;
%% Nfig: output figure number;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ccode = [0.5,0.5,0.5;0,0,1;0,1,0;1,0,0;0,0,0];

figure(Nfig)
maximize(Nfig)
    colormap(ccode);
    scatter(DNA_profile,RNA_profile(:,2),32,RNA_profile(:,1)+1,'fill')
    set(gca,'CLim',[1,5])
    xlabel('DNA intensity (A.U.)')
    ylabel('TX level (#)')
    hc = colorbar;
    set(hc,'YTick',[1:5],'YTickLabel',[0:4])
    title([image_folder,': cycle ',num2str(N_cycle)],'Interpreter','none')