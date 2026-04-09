%% live data change to tif; output folder: result'index', such as result1
bhh_path='X:/baohuihan/FISHIF-Time/Confocal/';
[num,txt,raw]=xlsread([bhh_path,'Livedata.xlsx'],'cycle12');%% input Embryo information in Excel
for j=1:length(txt(:,1))
    mkdir([bhh_path,char(txt(j,1)),'/result',sprintf('%d',num(j,4))])
img=tiffread([bhh_path,char(txt(j,1)),'/',char(txt(j,2)),'.lsm']);
Tstart=num(j,6);
Tend=num(j,6)+num(j,1)-1;
zlayer=num(j,5);
for time=Tstart:Tend
    for z=1:zlayer
        a=(time-1)*zlayer+z;
        Img=img(a).data;
        imwrite(Img,[bhh_path,char(txt(j,1)),'/result',sprintf('%d',num(j,4)),'/result',sprintf('%d',num(j,4)),'_',sprintf('%03d',z),'_',sprintf('%03d',time),'.tiff'],'tiff');
    end
end
end