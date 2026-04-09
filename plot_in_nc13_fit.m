%% spatiotemporal Bcd Hb Kr and fit in nc13 (fit parameter result:co, saved in "fit")
clear 
close all
[num,txt,raw]=xlsread('X:/baohuihan/FISHIF-Time/Confocal/Fit.xlsx','cycle13');%% input:data information integrated in Excel
protein_add1 = 'Bcd protein';
RNA_add = 'Kr RNA';
for j=2:length(txt(:,1))
    if char(txt(j,5))=='H'
       protein_add2 = 'Hb protein';
    else
        protein_add2 = 'Gt protein';
    end
bhh_path=char(txt(j,1));
    if ismissing(txt(j,2))== 1
 load([bhh_path,sprintf('%d',num(j,1)),'/fit/',char(txt(j,3)),'_',RNA_add,'_intensity_new3.mat']);
 load([bhh_path,sprintf('%d',num(j,1)),'/fit/',char(txt(j,3)),'_',protein_add1,'_intensity_new3.mat']);
 load([bhh_path,sprintf('%d',num(j,1)),'/fit/',char(txt(j,3)),'_',protein_add2,'_intensity_new3.mat']);
    else
       load([bhh_path,char(txt(j,2)),'/fit/',char(txt(j,3)),'_',RNA_add,'_intensity_new3.mat']);
       load([bhh_path,char(txt(j,2)),'/fit/',char(txt(j,3)),'_',protein_add1,'_intensity_new3.mat']);
       load([bhh_path,char(txt(j,2)),'/fit/',char(txt(j,3)),'_',protein_add2,'_intensity_new3.mat']);
    end
    T(j-1)=num(j,6);%% inferred time
    R(:,j-1)=fiout(:,2);
    P1(:,j-1)=nu_out(:,2);%Hb Gt
    P2(:,j-1)=nu_out2(:,2);%Bcd
    Tr(:,j-1)=T(j-1)*ones(21,1);
    Tp(:,j-1)=T(j-1)*ones(21,1);
end
PHb=P1;
PGt=P1;
THb=Tp;
TGt=Tp;
for j=2:length(txt(:,1))
    if char(txt(j,5))=='H'
       PGt(:,j-1)=0;
       TGt(:,j-1)=0;
    else
        PHb(:,j-1)=0;
       THb(:,j-1)=0;
    end
    
end


    nucleus_bin = 1:0.5:13;
    average_radius = 0.5;
 bin_max = min(nucleus_bin+average_radius,13);
 bin_min = max(nucleus_bin-average_radius,0.1);
 dq=jet(14);
%figure;
 for I_bin = 1:length(nucleus_bin)
     %% RNA
        RR=R.*double((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin)));
        RR(find(isnan(RR)==1))=0;
        RR(:,all(RR==0,1))=[];
        fi0(:,I_bin) = mean(RR,2);        
        timer(I_bin)=mean(Tr((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
 
        fi1(:,I_bin)=std(RR,0,2);
        [max_active(I_bin),Rnu]=max(fi0(:,I_bin));
        std_RNA(I_bin)=fi1(Rnu,I_bin);
  
        
      %% Protein
        PPHb=PHb.*double((Tp >= bin_min(I_bin))&(Tp <= bin_max(I_bin)));
        PPHb(find(isnan(PPHb)==1))=0;
        PPHb(:,all(PPHb==0,1))=[];
        nuHb(:,I_bin) = mean(PPHb,2);
        timepHb(I_bin)=mean(THb((THb >= bin_min(I_bin))&(THb <= bin_max(I_bin))));
        
        PP2=P2.*double((Tp >= bin_min(I_bin))&(Tp <= bin_max(I_bin)));
        PP2(find(isnan(PP2)==1))=0;
        PP2(:,all(PP2==0,1))=[];
        nu2(:,I_bin) = mean(PP2,2);
        nustd2(:,I_bin) = std(PP2,0,2);  
        timep(I_bin)=mean(Tp((Tp >= bin_min(I_bin))&(Tp <= bin_max(I_bin))));
        
 end

 %% fit
 EL=0.2:0.05:0.7;
% [R1,a,b]=unique(fi0(5:15,:)','rows','stable');
% R1=R1';
% [Bcd,~,~]=unique(nu2(5:15,:)','rows','stable');
% Bcd=Bcd';
% [Hb,~,~]=unique(nuHb(5:15,:)','rows','stable');
% Hb=Hb';
% t=unique(timer,'stable');
R1=fi0(5:15,:);
Bcd=nu2(5:15,:);
Hb=nuHb(5:15,:);
t=timer;
save([bhh_path,'fit/Bcd_13.mat'],'Bcd');
save([bhh_path,'fit/Hb_13.mat'],'Hb');
save([bhh_path,'fit/R1_13.mat'],'R1');
save([bhh_path,'fit/t_13.mat'],'t');
save([bhh_path,'fit/Bcd_13_all.mat'],'nu2');
save([bhh_path,'fit/Hb_13_all.mat'],'nuHb');
save([bhh_path,'fit/R1_13_all.mat'],'fi0');
%%
Kr=R1;
%Kr=R1;
co0=[0 20 2 8 5 7 7 2 8 5 7 20 2 8 5 7 1 t(1) 5 8 11];
lb=[0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0.3 3 4 7 10];
ub=[10 35 20 20 10 10 20 20 20 10 10 30 20 20 10 10 5 6 6 10 13];
[co,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@myfun_cycle13,co0,EL,Kr,lb,ub);
conf=nlparci(co,residual,'jacobian',jacobian,"Alpha",0.9);%90% confidence
save([bhh_path,'fit/co_13.mat'],'co');
%%   Å·À­·¨


% t=t(2:8);
% R1=R1(:,2:8);
% Bcd=Bcd(:,2:8);
% Hb=Hb(:,2:8);
x=zeros(length(EL),length(t));
y=zeros(length(EL),length(t));
R=zeros(length(EL),length(t));
%co=[0 2 3 3 3 3 2];
for n=1:length(EL)
    R(n,1)=R1(n,1);
for i=1:length(t)-1
    x(n,i)=Bcd(n,i);
    y(n,i)=Hb(n,i);
    dR1=co(1)+(co(2).*((x(n,i).^co(5))./(co(3).^co(5)+x(n,i).^co(5))).*(co(4).^co(6))./(co(4).^co(6)+y(n,i).^co(6)))-co(17).*R(n,i);
    dR2=co(1)+(co(7).*((x(n,i).^co(10))./(co(8).^co(10)+x(n,i).^co(10))).*(co(9).^co(11))./(co(9).^co(11)+y(n,i).^co(11)))-co(17).*R(n,i);
    dR3=co(1)+(co(12).*((x(n,i).^co(15))./(co(13).^co(15)+x(n,i).^co(15))).*(co(14).^co(16))./(co(14).^co(16)+y(n,i).^co(16)))-co(17).*R(n,i);
    dR4=-co(17).*R(n,i);
    if t(i+1)<=co(18)
        R(n,i+1)=R(n,i);
    elseif t(i)<=co(18)&&t(i+1)>co(18)
        R(n,i+1)=R(n,i)+(t(i+1)-co(18)).*dR1;
    elseif t(i)>co(18)&&t(i+1)<=co(19)
       R(n,i+1)=R(n,i)+(t(i+1)-t(i)).*dR1;
    elseif t(i)<=co(19)&&t(i+1)>co(19)
        R(n,i+1)=R(n,i)+(co(19)-t(i)).*dR1+(t(i+1)-co(19)).*dR2;
    elseif t(i)>co(19)&&t(i+1)<=co(20)
       R(n,i+1)=R(n,i)+(t(i+1)-t(i)).*dR2;
    elseif t(i)<=co(20)&&t(i+1)>co(20)
        R(n,i+1)=R(n,i)+(co(20)-t(i)).*dR2+(t(i+1)-co(20)).*dR3;
    elseif t(i)>co(20)&&t(i+1)<=co(21)
       R(n,i+1)=R(n,i)+(t(i+1)-t(i)).*dR3;
    elseif t(i)<=co(21)&&t(i+1)>co(21)
        R(n,i+1)=R(n,i)+(co(21)-t(i)).*dR3+(t(i+1)-co(21)).*dR4;
    elseif t(i)>co(21)
       R(n,i+1)=R(n,i)+(t(i+1)-t(i)).*dR4;
    end
end
end
figure
for i=1:length(t)
    %figure
    subplot(5,5,i)
    scatter3(Bcd(:,i),Hb(:,i),R(:,i),'filled','MarkerFaceColor','r');
    hold on
    scatter3(Bcd(:,i),Hb(:,i),Kr(:,i),'filled','MarkerFaceColor','b');
    xlabel('Bcd');ylabel('hb');
    zlim([0 20])
    title(['T=',num2str(t(i))]);
end

%% surf 
Lnwindow = 0.25;
nu_L1 = 0.5:0.5:8.5;
nu_L2 = 0.25:0.25:6.5;
r1=[Bcd(:,1),Hb(:,1),R1(:,1)];
rkr=zeros(length(nu_L1),length(nu_L2));
 for Lcenter1 = 1:length(nu_L1)
     for Lcenter2 = 1:length(nu_L2)
    nu_map = (r1(:,1) >= nu_L1(Lcenter1)-Lnwindow) & (r1(:,1) <= nu_L1(Lcenter1)+Lnwindow) & (r1(:,2) >= nu_L2(Lcenter2)-Lnwindow) & (r1(:,2) <= nu_L2(Lcenter2)+Lnwindow);
    rkr(Lcenter1,Lcenter2) = mean(r1(nu_map,3));
     end
 end
 rkr(find(isnan(rkr)==1))=0;
 [hb,bcd]=meshgrid(nu_L2,nu_L1);
 figure
 mesh(bcd,hb,rkr);
 
 ela=4;
 elp=10;
 for i=1:length(t)-1
     
    for Lcenter1 = 1:length(nu_L1)
     for Lcenter2 = 1:length(nu_L2)
    
    dR1=co(1)+(co(2).*((nu_L1(Lcenter1).^co(5))./(co(3).^co(5)+nu_L1(Lcenter1).^co(5))).*(co(4).^co(6))./(co(4).^co(6)+nu_L2(Lcenter2).^co(6)))-co(17).*rkr(Lcenter1,Lcenter2);
    dR2=co(1)+(co(7).*((nu_L1(Lcenter1).^co(10))./(co(8).^co(10)+nu_L1(Lcenter1).^co(10))).*(co(9).^co(11))./(co(9).^co(11)+nu_L2(Lcenter2).^co(11)))-co(17).*rkr(Lcenter1,Lcenter2);
    dR3=co(1)+(co(12).*((nu_L1(Lcenter1).^co(15))./(co(13).^co(15)+nu_L1(Lcenter1).^co(15))).*(co(14).^co(16))./(co(14).^co(16)+nu_L2(Lcenter2).^co(16)))-co(17).*rkr(Lcenter1,Lcenter2);
    dR4=-co(17).*rkr(Lcenter1,Lcenter2);
    if t(i+1)<=co(18)
        rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2);
    elseif t(i)<=co(18)&&t(i+1)>co(18)
        rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(t(i+1)-co(18)).*dR1;
    elseif t(i)>co(18)&&t(i+1)<=co(19)
       rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(t(i+1)-t(i)).*dR1;
    elseif t(i)<=co(19)&&t(i+1)>co(19)
        rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(co(19)-t(i)).*dR1+(t(i+1)-co(19)).*dR2;
    elseif t(i)>co(19)&&t(i+1)<=co(20)
       rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(t(i+1)-t(i)).*dR2;
    elseif t(i)<=co(20)&&t(i+1)>co(20)
        rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(co(20)-t(i)).*dR2+(t(i+1)-co(20)).*dR3;
    elseif t(i)>co(20)&&t(i+1)<=co(21)
       rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(t(i+1)-t(i)).*dR3;
    elseif t(i)<=co(21)&&t(i+1)>co(21)
        rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(co(21)-t(i)).*dR3+(t(i+1)-co(21)).*dR4;
    elseif t(i)>co(21)
       rkr(Lcenter1,Lcenter2)=rkr(Lcenter1,Lcenter2)+(t(i+1)-t(i)).*dR4;
    end
    
     end
    end
    pn1(:,i)=polyfit(Bcd(ela:elp,i),Bcd(ela:elp,i+1),1);
    pn2(:,i)=polyfit(Hb(ela:elp,i),Hb(ela:elp,i+1),1);
    nu_L1=pn1(1,i).*nu_L1+pn1(2,i);
    nu_L2=pn2(1,i).*nu_L2+pn2(2,i);
    [hb,bcd]=meshgrid(nu_L2,nu_L1);
    figure
%     mesh(bcd,hb,rkr,'FaceAlpha',0.2);
    surf(bcd,hb,rkr,'FaceAlpha',0.2);
    colormap(jet)
    hold on
%     scatter3(Bcd(ela:elp,i+1),Hb(ela:elp,i+1),R(ela:elp,i+1),'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4);
%     hold on
    scatter3(Bcd(ela:elp,i+1),Hb(ela:elp,i+1),R1(ela:elp,i+1),'filled','MarkerFaceColor','b','MarkerFaceAlpha',1);
    hold on
    for n=ela:elp
    plot3([Bcd(n,i+1),Bcd(n,i+1)],[Hb(n,i+1),Hb(n,i+1)],[R1(n,i+1),R(n,i+1)],'k');
    hold on
    end
    title(['T=',num2str(t(i+1))]);
    xlabel('Bcd');ylabel('hb');zlabel('kr');
    grid off
    
%     xlim([0 20])
%     ylim([0 20])
    zlim([0 25])
%     set(gca,'xtick',0:10:20);
%     set(gca,'ytick',0:10:20);
    set(gca,'ztick',0:10:20);
    view([-80 42])
 end