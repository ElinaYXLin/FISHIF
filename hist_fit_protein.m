function hist_fit_protein(lf_name)
%clear all
close all

sub_folder = 'Histogram_protein_A/';
sub_folder2 = 'Histogram_protein_M/';
sub_folder3 = 'Histogram_protein_P/';
% sub_folder_nullo = 'Histogram_default/';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
data_tail = '_raw.xls';
data_tail2 = '_raw.xlsx';
output_add = '_spot_fit_2';
mask_area_add = '_mask_area';
purify_add = '_purify';
mat_tail = '.mat';
fig_tail = '.fig';
hist_min = 0; 
hist_max = 2e5;%2e5
hist_bin = 2e3;%1e3
binw2 = 1e3;%0.5e3
% % fit_initial = [0.02,0.1,0.1,0.01,0.01,1e4,1e4,2e4,2e4,0];
% % fit_lower = [0,0,0,0,0,0,0,0,0,0];
% % fit_upper = [0.3,1,1,1,1,2e4,2e4,4e4,4e4,0];
% fit_initial = [0.02,0.1,0.1,0.01,0.01,0.01,0.8e4,0.8e4,1.6e4,1.6e4,0];
% fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
% fit_upper = [1,1,1,1,1,1,3e4,3e4,5e4,5e4,0];
% fit_initial = [0.01,0.04,0.04,0.001,0.001,0.001,1e4,4e3,2.5e4,2.5e4,0.01];
% fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
% fit_upper = [1,1,1,1,1,1,2e4,2e4,6e4,6e4,0.1];old
fit_initial = [0.01,0.04,0.04,0.001,0.001,0.001,2e4,2e4,0.6e5,0.6e5,0.01];
fit_lower = [0,0,0,0,0,0,0,0,0,0,0];
fit_upper = [1,1,1,1,1,1,5e4,5e4,2e5,2e5,0.1];
xlim0 = [0,2.5e4];%2e3
ylim0 = [0,5];
% Inten_lim = 500;

% bin3D = {0:5e1:1e5,0:0.1:10};
bin3D = {0:5e2:2.5e4,0:0.1:10};%4e3
rmax = 5;%6
relip = 0.95;%0.95

gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator
eval(['gau2 = @(x,xdata) ',gau2_gen(1)]);
fit_range = [10,10];
sigma0 = [1,500];
lb = [0,0,0,0,0,-inf,0];
ub = [10,100,bin3D{2}(end),10,bin3D{1}(end),inf,1];
options = optimset('Display','off');

is_purify = true;

%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name = lf_name;
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_hist = hist_min:hist_bin:hist_max;
x_hist_min = x_hist-binw2;
x_hist_max = x_hist+binw2;

for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    im_folder = folder_list{list_I,1};
        
    for list_J = 2%1:M1
        data_name = sub_list{list_J,3}(1:end-1);
        if exist([im_folder,sub_folder,data_name,data_tail],'file')
            raw_data = xlsread([im_folder,sub_folder,data_name,data_tail]);
        else
            raw_data = xlsread([im_folder,sub_folder,data_name,data_tail2]);
        end
        
        if exist([im_folder,sub_folder2,data_name,data_tail],'file')
            raw_data2 = xlsread([im_folder,sub_folder2,data_name,data_tail]);
        else
            raw_data2 = xlsread([im_folder,sub_folder2,data_name,data_tail2]);
        end
        
        if exist([im_folder,sub_folder3,data_name,data_tail],'file')
            raw_data3 = xlsread([im_folder,sub_folder3,data_name,data_tail]);
        else
            raw_data3 = xlsread([im_folder,sub_folder3,data_name,data_tail2]);
        end
        
        if is_purify && exist([im_folder,sub_folder,data_name,purify_add,mat_tail],'file')
            S = load([im_folder,sub_folder,data_name,purify_add,mat_tail]);
            I_purify = S.spot_true;
        else
            I_purify = true(size(raw_data,1),1);
        end
        
        if is_purify && exist([im_folder,sub_folder2,data_name,purify_add,mat_tail],'file')
            S = load([im_folder,sub_folder2,data_name,purify_add,mat_tail]);
            I_purify2 = S.spot_true;
        else
            I_purify2 = true(size(raw_data2,1),1);
        end
        
        if is_purify && exist([im_folder,sub_folder3,data_name,purify_add,mat_tail],'file')
            S = load([im_folder,sub_folder3,data_name,purify_add,mat_tail]);
            I_purify3 = S.spot_true;
        else
            I_purify3 = true(size(raw_data3,1),1);
        end
        
        raw_data = spfilter(raw_data(I_purify,:),rmax,rmax,relip,0);
        raw_data2 = spfilter(raw_data2(I_purify2,:),rmax,rmax,relip,0);
        raw_data3 = spfilter(raw_data3(I_purify3,:),rmax,rmax,relip,0);
        
        raw_data = cat(1,raw_data,raw_data2,raw_data3);
        X = [raw_data(:,1),sqrt(raw_data(:,2).*raw_data(:,3))];
        X2 = [raw_data2(:,1),sqrt(raw_data2(:,2).*raw_data2(:,3))];
        X3 = [raw_data3(:,1),sqrt(raw_data3(:,2).*raw_data3(:,3))];
% % %         X = [raw_data(:,1),raw_data(:,2)./raw_data(:,3)];
% % %         X2 = [raw_data2(:,1),raw_data2(:,2)./raw_data2(:,3)];
% % %         X3 = [raw_data3(:,1),raw_data3(:,2)./raw_data3(:,3)];
% % % %         load([im_folder,sub_folder,sub_list{list_J,3}(1:end-1),mask_area_add,mat_tail])
% % % %         mask_area1 = nnz(embryo_mask);
        
% % % %         load([im_folder,sub_folder2,sub_list{list_J,3}(1:end-1),mask_area_add,mat_tail])
% % % %         mask_area2 = nnz(embryo_mask);
        
        figure
        maximize(gcf)
        a(1) = subplot(2,3,1);
%             hist3(X,'CdataMode','auto','EdgeColor','none','Ctrs',{0:5e1:1e5,0:0.1:10})
            [N00,C0] = hist3(X,'CdataMode','auto','EdgeColor','none','Ctrs',bin3D);
            imagesc(C0{1},C0{2},log10(N00)');
            axis xy
            colorbar
%             view(2)
            hold on
%             plot3(Inten_lim*ones(size(ylim0)),ylim0,1000*ones(size(ylim0)),'w-','LineWidth',2)
            xlim(xlim0)
            ylim(ylim0)
            xlabel('Peak intensity (A.U.)')
            ylabel('/sigma (pixel)')
            title('Anterior spots');
        a(2) = subplot(2,3,2);
%             hist3(X2,'CdataMode','auto','EdgeColor','none','Ctrs',{0:5e1:1e5,0:0.1:10})
            [N0,C0] = hist3(X2,'CdataMode','auto','EdgeColor','none','Ctrs',bin3D);
            imagesc(C0{1},C0{2},log10(N0)');
            axis xy
            colorbar
%             view(2)
            hold on
%             plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2)
% %             plot3(Inten_lim*ones(size(ylim0)),ylim0,1000*ones(size(ylim0)),'w-','LineWidth',2)
            xlim(xlim0)
            ylim(ylim0)
            xlabel('Peak intensity (A.U.)')
            ylabel('/sigma (pixel)')
            title('Median spots');
        a(3) = subplot(2,3,3);
%             hist3(X2,'CdataMode','auto','EdgeColor','none','Ctrs',{0:5e1:1e5,0:0.1:10})
            [N0,C0] = hist3(X3,'CdataMode','auto','EdgeColor','none','Ctrs',bin3D);
            imagesc(C0{1},C0{2},log10(N0)');
            axis xy
            colorbar
%             view(2)
            hold on
%             plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2)
% %             plot3(Inten_lim*ones(size(ylim0)),ylim0,1000*ones(size(ylim0)),'w-','LineWidth',2)
            xlim(xlim0)
            ylim(ylim0)
            xlabel('Peak intensity (A.U.)')
            ylabel('/sigma (pixel)')
            title('Posterior spots');
        linkaxes(a)

        t_gate0 = false;
        h1 = [];
        while ~t_gate0
            if ~exist('Inten_lim') || isempty(Inten_lim)
                xlim(xlim0)
                ylim(ylim0)
                xlabel('Peak intensity (A.U.)')
                ylabel('/sigma (pixel)')
                title('Please pick up the boundary');
                [x0,y0] = ginput(2);
                if length(x0) > 1
                    p0 = polyfit(y0,x0,1); 
                else
                    p0(1) = 0; p0(2) = x0;
                end
            else
                p0(1) = 0; p0(2) = Inten_lim;
                t_gate0 = true;
            end

            delete(h1);
            axes(a(1))
                hold on
                h1(1) = plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2);
    %             xlim(xlim0)
            axes(a(2))
                hold on
                h1(2) = plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2);
    %             xlim(xlim0)
            axes(a(3))
                hold on
                h1(3) = plot(polyval(p0,ylim0),ylim0,'w-','LineWidth',2);
    %             xlim(xlim0)
            
            if ~t_gate0
                answer0 = questdlg('Is the gating OK','Gating confirmation','Yes','No','Yes');
                t_gate0 = strcmpi('Yes',answer0);
            end
        end
        
        [xxx,yyy] = meshgrid(C0{2},C0{1});
        Itrue0 = (polyval(p0,xxx(:)) <= yyy(:)) & (xxx(:) >= ylim0(1)) & (xxx(:) <= ylim0(2)) & (yyy(:) >= xlim0(1)) & (yyy(:) <= xlim0(2));
        [~,I0] = max(N00(Itrue0));
        Ntrue0 = find(Itrue0);
        rmax0 = xxx(Ntrue0(I0));
        Imax0 = yyy(Ntrue0(I0));
        bmax = Imax0*rmax0^2*2*pi;
        
        [I0i,I0j] = ind2sub(size(xxx),Ntrue0(I0));
        Itrue1 = false(size(xxx)); 
        Itrue1(max(I0i-fit_range(1),1):min(I0i+fit_range(1),size(Itrue1,1)),max(I0j-fit_range(2),1):min(I0j+fit_range(2),size(Itrue1,2))) = true;
        Itrue2 = Itrue1(:) & Itrue0;
        xxx1 = xxx(Itrue2);
        yyy1 = yyy(Itrue2);
        N11 = N00(Itrue2);
        xdata = [xxx1(:),yyy1(:)];
        ydata = N11(:); ydata = ydata/max(ydata);
        xf0 = [1, 1./2./sigma0(1).^2, rmax0, 1./2./sigma0(2).^2, Imax0, 0, 0];
        
        [xf,~,~,exitflag,~] = lsqcurvefit(gau2,xf0,xdata,ydata,lb,ub,options);
        
        rmax2 = xf(3);
        Imax2 = xf(5);
        % if exitflag > 0
        %     rmax2 = xf(3);
        %     Imax2 = xf(5);
        % else
        %     rmax2 = rmax;
        %     Imax2 = Imax;
        % end
        bmax2 = Imax2*rmax2^2*2*pi;

%% Fit histogram: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Itrue = polyval(p0,X(:,2)) <= X(:,1);
        Inten_spot = prod(raw_data(Itrue,1:3),2)*2*pi;
        n_hist = hist0(Inten_spot,x_hist_min,x_hist_max);

% % %         Itrue = polyval(p0,X2(:,2)) <= X2(:,1);
% % %         Inten_spot = prod(raw_data2(Itrue,1:3),2)*2*pi;
% % %         n_hist = hist0(Inten_spot,x_hist_min,x_hist_max);
        
% % %         Itrue2 = raw_data2(:,1) >= Inten_lim;
% % %         Inten_spot2 = prod(raw_data2(Itrue2,1:3),2)*2*pi;
% % %         n_hist2 = hist0(Inten_spot2,x_hist_min,x_hist_max);
% % %         
% % %         Itrue3 = raw_data3(:,1) >= Inten_lim;
% % %         Inten_spot3 = prod(raw_data3(Itrue3,1:3),2)*2*pi;
% % %         n_hist3 = hist0(Inten_spot3,x_hist_min,x_hist_max);
        
%         n_hist = conv(n_hist,[0.5,1,0.5]/2,'same');
% %         x_hist = x_hist(1:end-1);
% %         y_hist = n_hist(1:end-1)/sum(n_hist(1:end-1));
        y_hist = n_hist/sum(n_hist);%-n_hist2/sum(n_hist2)/4;
        mgau = @(a0,a1,a2,a3,a4,a5,b0,c0,b,c,d,x) a0*exp(-(x-b0).^2./2./c0.^2)+a1*exp(-(x-b).^2./2./c.^2)+a2*exp(-(x-2*b).^2./2./2./c.^2)+a3*exp(-(x-3*b).^2./3./2./c.^2)+a4*exp(-(x-4*b).^2./4./2./c.^2)+a5*exp(-(x-5*b).^2./5./2./c.^2)+d;
        spot_fit = fit( x_hist',y_hist',mgau,'StartPoint',fit_initial,'Lower',fit_lower,'Upper',fit_upper);
        b = spot_fit.b;
        save([im_folder,sub_folder,data_name,output_add,mat_tail],'spot_fit','x_hist','y_hist','b','bmax','bmax2');
        
        subplot(2,3,4:6)
            bar(x_hist,y_hist,'b')
            hold on
            x_fit = [hist_min:(hist_max-hist_min)/1000:hist_max];
%             y_fit = mgau(spot_fit.a0,spot_fit.a1,spot_fit.a2,spot_fit.a3,spot_fit.a4,spot_fit.b0,spot_fit.c0,spot_fit.b,spot_fit.c,spot_fit.d,x_fit);
            y_fit = feval(spot_fit,x_fit);
            plot(x_fit,y_fit,'r')
            xlabel('Spot intensity (A.U.)')
            ylabel('Frequency')
            title(['Spot intensity histogram fit: ',data_name,', b = ',num2str(b),', bmax = ',num2str(bmax),', bmax2 = ',num2str(bmax2)],'Interpreter','none')
        saveas(gcf,[im_folder,sub_folder,data_name,output_add,fig_tail])
        
%         b = 0;
%         if exist([im_folder,sub_folder_nullo]) ~= 7
%             mkdir([im_folder,sub_folder_nullo]);
%         end
%         save([im_folder,sub_folder_nullo,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
    end
end