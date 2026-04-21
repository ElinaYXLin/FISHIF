clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to manually remove bad nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_name = '02122013\';
result_folder = 'Results\';
mask_folder = 'masks\';
mask_name = 'mask.mat';
qmask_name = 'quick_mask2.mat';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading and plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = dir([folder_name,result_folder,'*.mat']);
for ii = 4:length(fname)
    data_name = fname(ii).name;
%     load([folder_name,result_folder,data_name],'nucleus_protein_profile','nucleus_protein2_profile')
    load([folder_name,result_folder,data_name],'nucleus_protein_profile','nucleus_GFP_profile')
    nucleus_protein2_profile = nucleus_GFP_profile;
    
    load([folder_name,mask_folder,data_name(1:end-4),'\',mask_name],'mask_stack')
    n1 = max(nucleus_protein_profile(:,2));
    v1 = max(nucleus_protein_profile(:,3));
    n2 = max(nucleus_protein2_profile(:,2));
    v2 = max(nucleus_protein2_profile(:,3));
    
    %%% plot protein1 profile
    figure(1)
    ha(1) = subaxis(2,2,1,1, 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(nucleus_protein_profile(:,1),nucleus_protein_profile(:,2),'b.')
        hold on
        xlabel('EL')
        ylabel('Protein 1')
        title([folder_name,result_folder,data_name])
    ha(2) = subaxis(2,2,1,2, 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(nucleus_protein_profile(:,2),nucleus_protein_profile(:,3),'b.')
        hold on
        xlabel('Mean 1')
        ylabel('Variance 1')
        title([folder_name,result_folder,data_name])
    ha(3) = subaxis(2,2,2,1, 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(nucleus_protein2_profile(:,1),nucleus_protein2_profile(:,2),'b.')
        hold on
        xlabel('EL')
        ylabel('Protein 2')
        title([folder_name,result_folder,data_name])
    ha(4) = subaxis(2,2,2,2, 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(nucleus_protein2_profile(:,2),nucleus_protein2_profile(:,3),'b.')
        hold on
        xlabel('Mean 1')
        ylabel('Variance 2')
        title([folder_name,result_folder,data_name])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Bad data point recognition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tn = true;   %%% switch1 to run "Bad data point recognition"
    tp = true;   %%% switch2 to run "Bad data point recognition"
    I_bad = zeros(0);
    while tn
        while tp
            I_temp = zeros(0);
            [x,y,n] = ginput(1);

            if ismember(n,[1:3])
                if gca == ha(1)
                    [~,I_temp] = min(pdist2([x,y/n1],[nucleus_protein_profile(:,1),nucleus_protein_profile(:,2)/n1]));
                elseif gca == ha(2)
                    [~,I_temp] = min(pdist2([x/n1,y/v1],[nucleus_protein_profile(:,2)/n1,nucleus_protein_profile(:,3)/v1]));
                elseif gca == ha(3)
                    [~,I_temp] = min(pdist2([x,y/n2],[nucleus_protein2_profile(:,1),nucleus_protein2_profile(:,2)/n2]));
                elseif gca == ha(4)
                    [~,I_temp] = min(pdist2([x/n2,y/v2],[nucleus_protein2_profile(:,2)/n2,nucleus_protein2_profile(:,3)/v2]));
                end

                if ~isempty(I_temp) && ~ismember(I_temp,I_bad)
                    I_bad = union(I_bad,I_temp);
                elseif ~isempty(I_temp) && ismember(I_temp,I_bad)
                    I_bad = setdiff(I_bad,I_temp);
                end

                axes(ha(1)); ho = get(ha(1),'Children'); if numel(ho) > 1; delete(ho(1)); end
                    plot(nucleus_protein_profile(I_bad,1),nucleus_protein_profile(I_bad,2),'r.')
                axes(ha(2)); ho = get(ha(2),'Children'); if numel(ho) > 1; delete(ho(1)); end
                    plot(nucleus_protein_profile(I_bad,2),nucleus_protein_profile(I_bad,3),'r.')
                axes(ha(3)); ho = get(ha(3),'Children'); if numel(ho) > 1; delete(ho(1)); end
                    plot(nucleus_protein2_profile(I_bad,1),nucleus_protein2_profile(I_bad,2),'r.')
                axes(ha(4)); ho = get(ha(4),'Children'); if numel(ho) > 1; delete(ho(1)); end
                    plot(nucleus_protein2_profile(I_bad,2),nucleus_protein2_profile(I_bad,3),'r.')
            else
                tp = false;
            end
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Display bad nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_bad = unique(I_bad);
        lmask = logical(mask_stack);
        bmask = ismember(mask_stack,I_bad);
        b2D = any(bmask,3);
        gmask = lmask & (~bmask);
        g2D = any(gmask,3); clear gmask bmask lmask;

        im0(:,:,1) = double(b2D);
        im0(:,:,3) = double(g2D);
        figure(2)
            imshow(im0); clear im0
            title([folder_name,result_folder,data_name,': y/n?'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Resave nuclei mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,~,n] = ginput(1);

        if n ~= 27
            Inu0 = 1:max(mask_stack(:));
            Inu1 = Inu0;
            Inu1(I_bad) = 0;
            Inu1(~ismember(Inu0,I_bad)) = 1:(length(Inu0)-length(I_bad));
            Inu1 = [0,Inu1];
            mask_stack = Inu1(mask_stack+1);
            bw_applied3D_save = logical(mask_stack);

            save([folder_name,mask_folder,data_name(1:end-4),'\',mask_name],'mask_stack','-append','-v7.3')
            save([folder_name,mask_folder,data_name(1:end-4),'\',qmask_name],'bw_applied3D_save','-v7.3')

            clear mask_stack bw_applied3D_save 
            close(1)
            close(2)
            
            tn = false;
        end
    end
end        
        

    






