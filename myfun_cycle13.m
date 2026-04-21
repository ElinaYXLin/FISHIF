%%   Å·À­·¨
function R=myfun_cycle13(co,EL)
load('Z:\kr-enhancer\fit2\t_13.mat');
load('Z:\kr-enhancer\fit2\Bcd_13.mat');
load('Z:\kr-enhancer\fit2\Hb_13.mat');
load('Z:\kr-enhancer\fit2\R1_13.mat');
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
end