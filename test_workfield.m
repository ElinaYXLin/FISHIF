Centroid_nucleus=Centroid_nucleus';
L=size(Centroid_nucleus,2);
Centroid_nucleus_chaos=Centroid_nucleus;
R=angle2dcm(pi/4,pi/3,pi/6); 
Centroid_nucleus_chaos=R*Centroid_nucleus_chaos;
Centroid_nucleus_chaos=Centroid_nucleus_chaos+[100 -200 3]';
Centroid_nucleus_chaos=Centroid_nucleus_chaos+[0*rand(2,L);0*rand(1,L)];


Centroid_nucleus=Centroid_nucleus(:,rand(1,L)>0.1);
Centroid_nucleus_chaos=Centroid_nucleus_chaos(:,rand(1,L)>0.1);

figure;
scatter3(Centroid_nucleus(1,:),Centroid_nucleus(2,:),Centroid_nucleus(3,:))
hold on
scatter3(Centroid_nucleus_chaos(1,:),Centroid_nucleus_chaos(2,:),Centroid_nucleus_chaos(3,:))

SignalAlignment3D(Centroid_nucleus',Centroid_nucleus_chaos')

figure;
scatter(Centroid_layers_focus(:,1),Centroid_layers_focus(:,2))
figure;
scatter(Centroid_nucleus(:,2),Centroid_nucleus(:,1))
