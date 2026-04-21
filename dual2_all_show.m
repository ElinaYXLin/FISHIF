function dual2_all_show(total_name,good_data,RNAI_out,RNAN_out,null_out,pro_center0,pro2_center0,pro_center,pro2_center,varargin)

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin)
    p_name = varargin{1}(1);
    p2_name = varargin{1}(2);
else
    p_name = 'Protein 1';
    p2_name = 'Protein 2';
end

if length(varargin) > 1
    fn = varargin{2};
else
    fn = [412,415,417];
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Mean plot of I regulation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RNAI_out0 = RNAI_out(:,:,good_data);
RNAI_value = RNAI_out0;
RNAI_value(isnan(RNAI_out0)) = 0;
RNAI_mean = sum(RNAI_value,3)./sum(~isnan(RNAI_out0),3);
RNAI_std0 = sqrt((sum(RNAI_value.^2,3)./sum(~isnan(RNAI_out0),3)-RNAI_mean.^2)./sum(~isnan(RNAI_out0),3));

figure(fn(1))
    clf
    subplot(1,2,1)
    imagesc(pro2_center0,pro_center0,RNAI_mean)
    axis xy
    title([total_name,': nucleus dual regulation (RNA level) curve (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['RNA level (#)'])

    subplot(1,2,2)
    imagesc(pro2_center0,pro_center0,RNAI_std0)
    axis xy
    title([total_name,': nucleus dual regulation (RNA level) error (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['RNA std (#)'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mean plot of N regulation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RNAN_out0 = RNAN_out(:,:,good_data);
RNAN_value = RNAN_out0;
RNAN_value(isnan(RNAN_out0)) = 0;
RNAN_mean = sum(RNAN_value,3)./sum(~isnan(RNAN_out0),3);
RNAN_std0 = sqrt((sum(RNAN_value.^2,3)./sum(~isnan(RNAN_out0),3)-RNAN_mean.^2)./sum(~isnan(RNAN_out0),3));

figure(fn(2))
    clf
    subplot(1,2,1)
    imagesc(pro2_center0,pro_center0,RNAN_mean)
    axis xy
    title([total_name,': nucleus dual regulation (foci #) curve (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['# of foci per nucleus'])

    subplot(1,2,2)
    imagesc(pro2_center0,pro_center0,RNAN_std0)
    axis xy
    title([total_name,': nucleus dual regulation (foci #) error (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['std of # of foci per nucleus'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mean plot of null rate: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
null_out0 = null_out(:,:,good_data);
null_value = null_out0;
null_value(isnan(null_out0)) = 0;
null_mean = sum(null_value,3)./sum(~isnan(null_out0),3);
null_std = sqrt(sum(null_value.^2,3)./sum(~isnan(null_out0),3)-null_mean.^2);
% dnull_mean = [diff(null_mean),nan]/(pro_center(2)-pro_center(1));
% dnull_std = [sqrt(null_std(1:end-1).^2+null_std(2:end).^2)/2,nan];
pshow_max = 5e-8;

figure(fn(3))
    clf
    subplot(1,2,1)
    imagesc(pro2_center,pro_center,null_mean)
    axis xy
    title([total_name,': foci active rate (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['Active rate'])

    subplot(1,2,2)
    imagesc(pro2_center,pro_center,null_std)
    axis xy
    title([total_name,': foci active rate std (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['Active rate std'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
