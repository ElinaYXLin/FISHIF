function Inten_compare(lf_name)
close all

% nfolder = '07022012/';
sub1 = 'Histogram_protein_A/';
sub2 = 'Histogram_protein_P/';
% fname = '06232012_2_OreR_FISHIF_hb48TMR_act5cA647_BcdSCA488_60X_2_001_new_raw';
data_tail = '_raw.xls';
in_tail = '.xls';
out_name = '_APcmp';
out_tail = '.fig';
nbin = [0:2e3:1e6];
nbin0 = [0:1e2:1e5];

%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name = lf_name;
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% AP comparison: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    nfolder = folder_list{list_I};
    sub0 = dir([nfolder,sub1,'*',data_tail]);
    M1 = length(sub0);
    data_name = {sub0.name};

    for list_J = 1:M1
        fname = data_name{list_J}(1:end-4);
        raw_A0 = xlsread([nfolder,sub1,fname,in_tail]);
        raw_P0 = xlsread([nfolder,sub2,fname,in_tail]);
        
%         ItrueA = sqrt(1-(raw_A(:,2)./raw_A(:,3)).^2) <= 0.9;
%         ItrueP = sqrt(1-(raw_P(:,2)./raw_P(:,3)).^2) <= 0.9;
%         ItrueA = raw_A0(:,8) >= 4 & raw_A0(:,8) <= 16 & raw_A0(:,1) >= 0;
%         ItrueP = raw_P0(:,8) >= 4 & raw_P0(:,8) <= 16 & raw_P0(:,1) >= 0;
%         raw_A = raw_A0(ItrueA,:);
%         raw_P = raw_P0(ItrueP,:);
        raw_A = raw_A0;
        raw_P = raw_P0;

        IntenA = raw_A(:,1).*raw_A(:,2).*raw_A(:,3)*2*pi;
        IntenP = raw_P(:,1).*raw_P(:,2).*raw_P(:,3)*2*pi;
        nIA = hist(IntenA,nbin);
        nIP = hist(IntenP,nbin);

        IntenA0 = raw_A0(:,1);
        IntenP0 = raw_P0(:,1);
        nIA0 = hist(IntenA0,nbin0);
        nIP0 = hist(IntenP0,nbin0);
        
        figure(1)
        clf
        subplot(2,1,1)
            plot(nbin,nIA,'r','DisplayName','Anterior')
            hold on
            plot(nbin,nIP,'b','DisplayName','Posterior')
            plot(nbin,nIA-nIP,'g','DisplayName','Difference')
            xlabel('Spot Intensity (A.U.)')
            ylabel('Count (#)')
            title([nfolder,fname],'Interpreter','none')
            legend('Show')
            xlim([0,1e5])
            
        subplot(2,1,2)
            plot(nbin0,nIA0,'r','DisplayName','Anterior')
            hold on
            plot(nbin0,nIP0,'b','DisplayName','Posterior')
            plot(nbin0,nIA0-nIP0,'g','DisplayName','Difference')
            xlabel('Spot Peak (A.U.)')
            ylabel('Count (#)')
            title([nfolder,fname],'Interpreter','none')
            legend('Show')
            xlim([0,1e4])
        
        saveas(1,[nfolder,sub1,fname,out_name,out_tail])
    end
end
