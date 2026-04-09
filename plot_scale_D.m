%% average live diameter trend
clear
bhh_path='X:/baohuihan/FISHIF-Time/Confocal/';
[num,txt,raw]=xlsread([bhh_path,'Livedata.xlsx'],'cycle12');%% input Embryo information in Excel
Tadd=0;
figure
for j=1:length(txt(:,1))
    if char(txt(j,5))=='F'
        continue
    else

load([bhh_path,char(txt(j,1)),'/Diameter',sprintf('%d',num(j,4)),'_cycle12new.mat']);
for i=1:length(Dout)
    %% for nc11
%     if num(j,11)==0.6
%     T(i,j)=(i/length(Dout)*7)-0.4;
%     else
%         T(i,j)=i/length(Dout)*7;
%     end
    %% for nc12-13 (Duration = 9 and 13)
    Duration=9;
    T(i,j)=i/length(Dout)*Duration;
    D(i,j)=Dout(i,1);    
end
% for i=1:length(Dout_th)
%     T(i,j)=i;
%     D(i,j)=Dout_th(i,1);    
% end
 scatter(T(:,j)+Tadd,D(:,j),'filled','MarkerFaceColor',[0.5 0.8 0.5],'MarkerFaceAlpha',.4);
 hold on
    end
 
end
nucleus_bin = 1:1:Duration;% bin 1
 average_radius = 0.5;% nc14:1
 bin_max = min(nucleus_bin+average_radius,Duration);
 bin_min = max(nucleus_bin-average_radius,0.1);
 fi0 = zeros(size(nucleus_bin));
 fi1 = zeros(size(nucleus_bin));
 for I_bin = 1:length(nucleus_bin)
        fi0(I_bin) = mean(D((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        fi1(I_bin) = std(D((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
 end
 nucleus_bin=nucleus_bin+Tadd;
 %figure
 errorbar(nucleus_bin,fi0,fi1,'b');
 %plot(nucleus_bin,fi0,'r')
 xlabel('T/min');
 ylabel('Diameter');
 title('23650 live embryo(meanD)');
 hold on
fimax=fi0+fi1;
fimin=fi0-fi1;
 for i = 1:length(nucleus_bin)-1
    x = [nucleus_bin(i),nucleus_bin(i+1),nucleus_bin(i+1),nucleus_bin(i)];
    y = [fimax(i),fimax(i+1),fimin(i+1),fimin(i)];
    h1 = fill(x,y,'m');
    set(h1,'Facecolor',[0.8 0.9 0.8],'FaceAlpha',0.3,'EdgeColor','none');
 end
 hold on
 
Tout=T(:);
Tout(Tout==0)=[];
Dout=D(:);
Dout(Dout==0)=[];
%% fixed diamter and infer scaling magnitude
 bhh_path='X:/baohuihan/FISHIF-Time/Confocal/';
 [num2,txt2,raw2]=xlsread([bhh_path,'Fixeddata.xlsx'],'cycle12');%% input fixed data information; with Diameter and inferred time
 hold on;plot(num2(:,5),num2(:,4)/0.83,'o')
fixT=num2(:,5)+Tadd;
fixD=num2(:,4);
% Traw = zeros(size(nucleus_bin));
% Draw = zeros(size(nucleus_bin));
% %  fi2 = zeros(size(nucleus_bin));
%  for I_bin = 1:length(nucleus_bin)
%         Traw(I_bin) = mean(fixT((fixT >= bin_min(I_bin))&(fixT <= bin_max(I_bin))));
%         Draw(I_bin) = mean(fixD((fixT >= bin_min(I_bin))&(fixT <= bin_max(I_bin))));
% %         fi2(I_bin) = std(fixD((fixT >= bin_min(I_bin))&(fixT <= bin_max(I_bin))));
%  end
% hold on; plot(Traw,Draw,'o')
 %% automatic infer 1
% Draw=rmmissing(Draw');
% Traw=rmmissing(Traw');
 Draw=rmmissing(num2(:,4));
 Traw=rmmissing(num2(:,5)+Tadd);
 Tout=repmat(Traw,1,length(nucleus_bin));
 Tbin=repmat(nucleus_bin,length(Traw),1);
 Tres=Tout-Tbin;
 rowsAllNegative = all(Tres < 0, 2);
 Tres(rowsAllNegative,1)=0;
 Tres(Tres<0)=100;
 [min_Res,index]=min(Tres,[],2);
 Dbin0=fi0;
 for n=1:length(index)
 Dbin(n)=((Dbin0(min(max(index),index(n)+1))-Dbin0(index(n)))/0.5*min_Res(n))+Dbin0(index(n));
 end
 for j=1:10
    rate(j)=0.71+0.02*j;
    L(j)=norm(Draw-Dbin'*rate(j));
end
Lout(:,1)=rate;
Lout(:,2)=L;
figure;plot(Lout(:,1),Lout(:,2));

%% manual check 2
for i=1:length(Traw)
    plot([Traw(i),Traw(i)],[22,55])
    hold on
end
[x,y]=ginput;% live
[x2,y2]=ginput;%fix

for j=1:10
    rate(j)=0.8+0.01*j;
    L(j)=norm(y2-y*rate(j));
end
Lout(:,1)=rate;
Lout(:,2)=L;
figure;plot(Lout(:,1),Lout(:,2));

 
 