
folder0 = '05102015/';
data_sub1 = 'Histogram/';
data_sub2 = 'Histogram_RNA2/';
data_tail = '*_raw.xls';

fname1 = dir([folder0,data_sub1,data_tail]);
fname2 = dir([folder0,data_sub2,data_tail]);

r = zeros(size(fname1));

for ii = 1:length(fname1)
    temp1 = xlsread([folder0,data_sub1,fname1(ii).name]);
    temp2 = xlsread([folder0,data_sub2,fname2(ii).name]); 
    
    Inten1 = prod(temp1(:,1:3),2);
    Inten2 = prod(temp2(:,1:3),2);
    
    r(ii) = mean(Inten1)/mean(Inten2);
end

r
    