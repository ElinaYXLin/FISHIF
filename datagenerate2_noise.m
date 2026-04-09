%% Background noise of DNA images and saturation out (output in DAPI;  "_nover" image output)
clear
bhh_path='X:/baohuihan/FISHIF-Time/Confocal/';

datagenerate_out(bhh_path,'A','s')
datagenerate_out(bhh_path,'M','s')
datagenerate_out(bhh_path,'P','s')
datagenerate_out(bhh_path,'','s')
datagenerate_out(bhh_path,'A','m')
datagenerate_out(bhh_path,'M','m')
datagenerate_out(bhh_path,'P','m')
datagenerate_out(bhh_path,'','m')
datagenerate_out(bhh_path,'A','l')
datagenerate_out(bhh_path,'M','l')
datagenerate_out(bhh_path,'P','l')
datagenerate_out(bhh_path,'','l')

function datagenerate_out(bhh_path,location,rect)%% input cropped mode: 's' 'm' 'l';input cropped mode: 'A' 'M' 'P' ''

[num,txt,raw]=xlsread([bhh_path,'Fixeddata.xlsx'],'cycle11');%% input Embryo information in Excel
for n=1:length(txt)
     bhh_path1=[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_',rect,'_',char(location),'w/',char(txt(n,2))];
     bhh_path2=[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_',rect,'_',char(location),'w_n/',char(txt(n,2))];
     bhh_path3=[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_',rect,'_',char(location),'w_nover/',char(txt(n,2))];
mkdir(bhh_path2);
mkdir(bhh_path3);
files = dir(fullfile(bhh_path1,'*.tiff'));
for j=1
    files1 = dir(fullfile(bhh_path2,'*.tiff'));
for i=1:size(files,1)
    
I=imread([bhh_path1,'/',files(i).name]);
Iz=double(I);
Iz1=sort(Iz(:),'descend');
count=hist(Iz1,unique(Iz1));
back=unique(Iz1);
[In,b]=max(count);
bn=back(b);


fit_initial=[200,100,10,5,200,bn,5];
 fit_higher=[300,200,50,50,max(count)*1.2,120,50];%60
 fit_lower=[0,0,0,0,0,0,0];
 ft = fittype(@(a1,b,c,d,a0,b0,c0,x) a0*exp(-(x-b0).^2./2./c0.^2)+a1*exp(-(x-b).^2./2./c.^2)+d );

 fitresult= fit( back, count', ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower);
 bn2=fitresult.b0;
 c2=fitresult.c0;
 
 
if c2<20
    var_gauss=sqrt(20.^2-c2.^2);
else
    var_gauss=0;
end

%% imnoise add
if bn2<60
   %% gausian
    m=60-bn2;
    %var_gauss=sqrt(max(1,(20.^2-c2.^2)));
    I_noise_add=m+var_gauss*randn(size(Iz));
    I_noise2=Iz+I_noise_add;%noise1 0.1 noise2 0.15
    
    Izn=round(sort(I_noise2(:),'descend'));
    countn=hist(Izn,unique(Izn));
    backn=unique(Izn);
    [~,b2]=max(countn(1:(length(countn)-1)));
    bn_noiseadd=backn(b2);
    
    fit_initial2=[100,100,60,5,300,60,20];
    fit_higher2=[300,200,50,50,max(countn)*1.2,120,60];%60
    fit_lower2=[0,0,0,0,0,0,0];
    ft2 = fittype(@(a1,b,c,d,a0,b0,c0,x) a0*exp(-(x-b0).^2./2./c0.^2)+a1*exp(-(x-b).^2./2./c.^2)+d );
    fitresult2= fit( backn, countn', ft2, 'StartPoint', fit_initial2, 'Upper', fit_higher2, 'Lower', fit_lower2);
    bn_add=fitresult2.b0;
    c2_add=fitresult2.c0;
 
    %I_noise2=I_noise2/max(max(I_noise2));
else
    I_noise2=Iz;
end
%% noise add out
Iout2=uint8(I_noise2);
 imwrite(Iout2,[bhh_path2,'/',files(i).name]);

%% saturation
Izn1=sort(Iout2(:),'descend');
x=find(Izn1>=255);
rate(i)=length(x)/length(Izn1);
Izn2=Izn1(1:length(Izn1)*0.02);
norm(i)=min(Izn2);
Iout2z=double(Iout2);
Ioutnover=Iout2z/double(norm(i));
imwrite(Ioutnover,[bhh_path3,'/',files(i).name]);
 

 end
end
Rate(n)=mean(rate);
Norm(n)=mean(norm);
BN(n,1)=median(bn2);%% get noise peak
BN(n,2)=median(c2);

end
end

 
 

