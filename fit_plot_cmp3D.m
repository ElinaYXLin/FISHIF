function fit_plot_cmp3D(xf,xyz,xy_radius,imstack0)

gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2+x(',num2str(n),',8)*(xdata(:,3)-x(',num2str(n),',9)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator

textfun = 'gau2 = @(x,xdata) ';
for n = 1:size(xf,1)
    textfun = [textfun,gau2_gen(n),'-x(',num2str(n),',7)/',num2str(size(xf,1)),'*(',num2str(size(xf,1)),'-1)','+'];
end
textfun = [textfun(1:(end-1)),';'];
eval(textfun);

xmin = xyz(1)-xy_radius;
xmax = xyz(1)+xy_radius;
ymin = xyz(2)-xy_radius;
ymax = xyz(2)+xy_radius;
zmin = xyz(3);
zmax = xyz(3);
bin0 = 0.1;
bin2 = 1;

[yg,xg,zg] = meshgrid(ymin:ymax,xmin:xmax,xyz(3));
zzz = double(imstack0(xmin:xmax,ymin:ymax));


xdata = [xg(:),yg(:),zg(:)];
ydata = zzz(:);

[yg1,xg1,zg1] = meshgrid(ymin:bin0:ymax,xmin:bin0:xmax,xyz(3));
zzz1 = zeros(size(xg1));
zzz1(:) = gau2(xf,[xg1(:),yg1(:),zg1(:)]);
[yg2,xg2,zg2] = meshgrid(ymin:bin2:ymax,xmin:bin2:xmax,xyz(3));
zzz2 = zeros(size(xg2));
zzz2(:) = gau2(xf,[xg2(:),yg2(:),zg2(:)]);
zzz0 = zeros(size(xg));
zzz0(:) = gau2(xf,[xg(:),yg(:),zg(:)]);

figure
    hs = stem3(yg(:),xg(:),zzz(:),'Color','r','LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',5,'ShowBaseLine','off');
    hsc = get(hs,'Children');
    z0 = get(hsc(1),'ZData');
    z0(2:3:end) = zzz0(:);
    set(hsc(1),'ZData',z0);
    
    hold on
%     mesh([ymin:bin2:ymax],[xmin:bin2:xmax],zzz2,'EdgeColor',[0,0,0]);
    meshx([ymin:bin2:ymax],[xmin:bin2:xmax],[ymin:bin0:ymax],[xmin:bin0:xmax],zzz1,[0,0,0]);

    h = surf([ymin:bin0:ymax],[xmin:bin0:xmax],zzz1);
    colormap([0,1,0])
    set(h,'edgecolor','none','LineStyle','none')
    camlight
    lighting phong;
    alpha(0.4)

    xlim([ymin,ymax]);
    ylim([xmin,xmax]);
%     axis equal
%     box on
    box off
    grid on
    xlabel('Pixel','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel('Pixel','FontName','Arial','FontSize',16,'FontWeight','bold')
    zlabel('Intensity (A.U.)','FontName','Arial','FontSize',16,'FontWeight','bold')
    set(gca,'FontName','Arial','FontSize',16,'FontWeight','normal')
