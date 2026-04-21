function outimage = corr_dist(inimage,channel_name,resolution)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to correct the distortion of the objective: %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global standard_data
st_list = fieldnames(standard_data);
imsize = size(inimage);
outimage = uint16(zeros(imsize));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Channel matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I_channel = 1:imsize(3)
    channel_match = find(strncmp(channel_name{I_channel},st_list,length(channel_name{I_channel})),1);
    if ~isempty(channel_match)
        eval(['current_standard = standard_data.',st_list{channel_match},';']);
        
        %%% Fitting parameter collection: %%% =============================
        fit2D = current_standard.fit;
        resolution0 = current_standard.resolution;
        stsize = current_standard.lattice;
        %%% ===============================================================
        
        %%% Image correction mask generation: %%% =========================
        X_start = (stsize(1)*resolution0-imsize(1)*resolution)/2/resolution0+1;
        Y_start = (stsize(2)*resolution0-imsize(2)*resolution)/2/resolution0+1;
        [X,Y] = meshgrid((X_start+[0:(imsize(1)-1)]*resolution/resolution0),(Y_start+[0:(imsize(2)-1)]*resolution/resolution0));
        XY = [reshape(X,imsize(1)*imsize(2),1),reshape(Y,imsize(1)*imsize(2),1)];
        Z1 = polyvaln(fit2D,double(XY));
        imcorr = reshape(Z1,imsize(1),imsize(2));
        %%% ===============================================================
    else
        imcorr = ones(imsize(1:2));
    end
    outimage(:,:,I_channel) = uint16(double(inimage(:,:,I_channel))./imcorr);
end