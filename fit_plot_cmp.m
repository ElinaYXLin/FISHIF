function fit_plot_cmp(xf,x,y,xy_radius,imstack0)

gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator

textfun = 'gau2 = @(x,xdata) ';
for n = 1:size(xf,1)
    textfun = [textfun,gau2_gen(n),'-x(',num2str(n),',7)/',num2str(size(xf,1)),'*(',num2str(size(xf,1)),'-1)','+'];
end
textfun = [textfun(1:(end-1)),';'];
eval(textfun);

xmin = x-xy_radius;
xmax = x+xy_radius;
ymin = y-xy_radius;
ymax = y+xy_radius;
bin0 = 0.1;
bin2 = 1;

[yg,xg] = meshgrid(ymin:ymax,xmin:xmax);
zg = double(imstack0(xmin:xmax,ymin:ymax));


xdata = [xg(:),yg(:)];
ydata = zg(:);

[yg1,xg1] = meshgrid(ymin:bin0:ymax,xmin:bin0:xmax);
zg1 = zeros(size(xg1));
zg1(:) = gau2(xf,[xg1(:),yg1(:)]);
[yg2,xg2] = meshgrid(ymin:bin2:ymax,xmin:bin2:xmax);
zg2 = zeros(size(xg2));
zg2(:) = gau2(xf,[xg2(:),yg2(:)]);
zg0 = zeros(size(xg));
zg0(:) = gau2(xf,[xg(:),yg(:)]);

figure
    hs = stem3(yg(:),xg(:),zg(:),'Color','r','LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',5,'ShowBaseLine','off');
    hsc = get(hs,'Children');
    z0 = get(hsc(1),'ZData');
    z0(2:3:end) = zg0(:);
    set(hsc(1),'ZData',z0);
    
    hold on
%     mesh([ymin:bin2:ymax],[xmin:bin2:xmax],zg2,'EdgeColor',[0,0,0]);
    meshx([ymin:bin2:ymax],[xmin:bin2:xmax],[ymin:bin0:ymax],[xmin:bin0:xmax],zg1,[0,0,0]);

    h = surf([ymin:bin0:ymax],[xmin:bin0:xmax],zg1);
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
