function Rna_distance(Rectify,Rna1_path,Rna2_path,save_path)
data_list1 = readmatrix(Rna1_path);
data_list2 = readmatrix(Rna2_path);
rna_local_1=data_list1(:,[6 7]);
rna_local_2=data_list2(:,[6 7]);
figure
scatter(rna_local_1(:,2),rna_local_1(:,1),12)
hold on
scatter(rna_local_2(:,2),rna_local_2(:,1),12)
xlim([0 Rectify.size(2)])
ylim([0 Rectify.size(1)])
saveas(gcf,[save_path,'rna-local','.fig'])
saveas(gcf,[save_path,'rna-local','.png'])
close(gcf)
Distance=@(x,y) sum((rna_local_1(x,:)-rna_local_2(y,:)).^2,2).^0.5;
[X,Y] = ndgrid(1:size(rna_local_1,1),1:size(rna_local_2,1));
Distance_matrix=Distance(X(:),Y(:));
dis_max=max(Distance_matrix);
dis_min=min(Distance_matrix);
bin=0.2;Edge=0:bin:dis_max;
Distrbution=histcounts(Distance_matrix,Edge);
%for n=1:1:dis_max/bin
    %Distrbution(:,n)=Distrbution(:,n)/(pi*((n*bin)^2-((n-1)*bin)^2));
%end
%DistrbutionSum=cumsum(Distrbution)./(Edge(1:end-1)+bin/2).^2;
figure
bar(Edge(1:end-1)+bin/2,Distrbution,'hist')
xlim([1 dis_min*100])
xlabel('DistanceSum')
ylabel('Number')
Roi=drawpoint;
distance_close=Roi.Position(1);
saveas(gcf,[save_path,'Distance-distrbution','.fig'])
saveas(gcf,[save_path,'Distance-distrbution','.png'])
close(gcf)
Distance_matrix=reshape(Distance_matrix,size(X));
[closeIndex1,closeIndex2]=find(Distance_matrix<=distance_close);
data_list1_close=data_list1(closeIndex1,:);
data_list2_close=data_list2(closeIndex2,:);
rna_local_1=data_list1_close(:,[6 7]);
rna_local_2=data_list2_close(:,[6 7]);
figure
scatter(rna_local_1(:,2),rna_local_1(:,1),12)
hold on
scatter(rna_local_2(:,2),rna_local_2(:,1),12)
xlim([0 Rectify.size(2)])
ylim([0 Rectify.size(1)])
saveas(gcf,[save_path,'rna-local-close','.fig'])
saveas(gcf,[save_path,'rna-local-close','.png'])
close(gcf)


[~,endIndex]=regexp(Rna1_path,'\');
Rna1_name=Rna1_path(endIndex(end)+1:end);
Rna1_name_close=[Rna1_name(1:end-4),'_close.xls'];
[~,endIndex]=regexp(Rna2_path,'\');
Rna2_name=Rna2_path(endIndex(end)+1:end);
Rna2_name_close=[Rna2_name(1:end-4),'_close.xls'];
writematrix(data_list1_close,[save_path,Rna1_name_close],'WriteMode','overwrite');
writematrix(data_list2_close,[save_path,Rna2_name_close],'WriteMode','overwrite');
end
