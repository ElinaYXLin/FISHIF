function imcorr = corr_mask(imsize,channel_name,resolution,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to correct the distortion of the objective: %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global standard_data standard_data_60X

if ~isempty(varargin) && ~isempty(varargin{1})
    standard_data0 = standard_data_60X;
else
    standard_data0 = standard_data;
end
    
if length(varargin) >= 2 && ~isempty(varargin{2})
    fit_model = varargin{2};
else
    fit_model = 0;
end
    
st_list = fieldnames(standard_data0);
imcorr = zeros(imsize);

pp = 1e-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Channel matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I_channel = 1:imsize(3)
    channel_match = find(strncmp(channel_name{I_channel},st_list,length(channel_name{I_channel})),1);
    if ~isempty(channel_match)
        eval(['current_standard = standard_data0.',st_list{channel_match},';']);
        
        if fit_model == 0
            %%% Fitting parameter collection: %%% =============================
            fit2D = current_standard.fit;
            resolution0 = current_standard.resolution;
            stsize = current_standard.lattice;
            %%% ===============================================================

            %%% Image correction mask generation: %%% =========================
            X_start = (stsize(1)*resolution0-imsize(1)*resolution)/2/resolution0+1;
            Y_start = (stsize(2)*resolution0-imsize(2)*resolution)/2/resolution0+1;
            [X,Y] = meshgrid((Y_start+[0:(imsize(2)-1)]*resolution/resolution0),(X_start+[0:(imsize(1)-1)]*resolution/resolution0));
            XY = [reshape(X,imsize(1)*imsize(2),1),reshape(Y,imsize(1)*imsize(2),1)];
            Z1 = polyvaln(fit2D,double(XY));
            imcorr(:,:,I_channel) = reshape(Z1,imsize(1),imsize(2));
            %%% ===============================================================
            
        elseif fit_model == 1
            %%% Fitting parameter collection: %%% =============================
            im00 = current_standard.image0;
            resolution0 = current_standard.resolution;
            stsize = current_standard.lattice;
            %%% ===============================================================

            %%% Image correction mask generation: %%% =========================
            X00 = ((1:stsize(1))-(1+stsize(1))/2)*resolution0;
            Y00 = ((1:stsize(2))-(1+stsize(2))/2)*resolution0;
            
            X11 = ((1:imsize(1))-(1+imsize(1))/2)*resolution;
            Y11 = ((1:imsize(2))-(1+imsize(2))/2)*resolution;
            
            imcorr(:,:,I_channel) = csaps({X00,Y00},im00,pp,{X11,Y11});
            %%% ===============================================================
            
        else
            imcorr(:,:,I_channel) = ones(imsize(1:2));
        end
    else
        imcorr(:,:,I_channel) = ones(imsize(1:2));
    end
end