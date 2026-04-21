%%   Å·À­·¨
function R=myfun_cycle11(co,EL)
load('Z:\kr-enhancer\fit2\t_11.mat');
load('Z:\kr-enhancer\fit2\Bcd_11.mat');
load('Z:\kr-enhancer\fit2\Hb_11.mat');
load('Z:\kr-enhancer\fit2\R1_11.mat');
% t=t(3:5);
% R1=R1(:,3:5);
% Bcd=Bcd(:,3:5);
% Hb=Hb(:,3:5);
x=zeros(length(EL),length(t));
y=zeros(length(EL),length(t));
R=zeros(length(EL),length(t));
%co=[0 2 3 3 3 3 2];
for n=1:length(EL)
    R(n,1)=R1(n,1);
for i=1:length(t)-1
    x(n,i)=Bcd(n,i);
    y(n,i)=Hb(n,i);
    dR1=co(1)+(co(2).*((x(n,i).^co(5))./(co(3).^co(5)+x(n,i).^co(5))).*(co(4).^co(6))./(co(4).^co(6)+y(n,i).^co(6)))-co(7).*R(n,i);
    dR2=-co(7).*R(n,i);
    if t(i+1)<=co(8)
        R(n,i+1)=R(n,i);
    elseif t(i)<=co(8)&&t(i+1)>co(8)
        R(n,i+1)=R(n,i)+(t(i+1)-co(8)).*dR1;
    elseif t(i)>co(8)&&t(i+1)<=co(9)
       R(n,i+1)=R(n,i)+(t(i+1)-t(i)).*dR1;
    elseif t(i)<=co(9)&&t(i+1)>co(9)
        R(n,i+1)=R(n,i)+(co(9)-t(i)).*dR1+(t(i+1)-co(9)).*dR2;
    elseif t(i)>co(9)
       R(n,i+1)=R(n,i)+(t(i+1)-t(i)).*dR2;
    end
end
end
end