function EL_info = get_EL(em_mask,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to extract EL information from the embryo mask %%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EL_info (output): EL information (x0,y0,x1,y1,EL_length);
%% em_mask (input): embryo mask;
%% varargin (input); {t_al: alternative method for determining the poles};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) && ~isempty(varargin{1})
    t_al = varargin{1};
else
    t_al = false;
end

em_mask = bwconvhull(em_mask);
[yall,xall] = find(bwperim(em_mask));
xy_all = [xall,yall];

if size(em_mask,2)/size(em_mask,1) > 1.4 || (~t_al)
    distXY = pdist2(xy_all,xy_all);
    axis_length = max(max(distXY));
    [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
    x0 = xy_all(extreme0,1);
    y0 = xy_all(extreme0,2);
    x1 = xy_all(extreme1,1);
    y1 = xy_all(extreme1,2);
    L2_extreme = (x1-x0)^2+(y1-y0)^2;
else
    I1 = xy_all(:,1) < size(em_mask,2)/2;
    I2 = xy_all(:,1) >= size(em_mask,2)/2;
    if std(xy_all(I1,2)) <= std(xy_all(I2,2))
        [x0,I0] = min(xy_all(:,1));
        y0 = xy_all(I0,2);
        x1 = 2*size(em_mask,2)-x0;
        y1 = y0;
        axis_length = x1-x0;
        L2_extreme = axis_length^2;
    else
        [x0,I0] = max(xy_all(:,1));
        y0 = xy_all(I0,2);
        x1 = 2*1-x0;
        y1 = y0;
        axis_length = x1-x0;
        L2_extreme = axis_length^2;
    end
end

EL_info = [x0,y0,x1,y1,L2_extreme];
