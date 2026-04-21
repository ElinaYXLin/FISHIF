function hist_resave_RNA(lf_name)
%clear all
close all

sub_folder = 'Histogram/';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
load_tail = '_value.mat';
data_tail = '_raw.xls';
foci_add = '_foci';
th_add = '_th';
fig_tail = '.fig';

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
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    im_folder = folder_list{list_I,1};
        
    for list_J = 1:M1
        imname = sub_list{list_J,3}(1:end-1);
        load_name = [imname,load_tail];
        load([im_folder,sub_folder,load_name])
%         stack_th0 = stack_th;
%         clear stack_th
        
        figure(1)
            clf
            if isempty(t0)
                t0 = nan(size(lim));
            end
            if isempty(dt0)
                dt0 = nan(size(lim));
            end
            [AX,H1,H2] = plotyy(lim,t0,lim,dt0);
            hold on
            plot(stack_th0*[1,1],ylim,'k--')
            title(['Threshold - # profile and the selected threshold value (',imname,') threshold = ',num2str(stack_th0)],'Interpreter','none')
            xlabel('Peak height threshold value (A.U.)')
            set(get(AX(1),'Ylabel'),'String','spot #') 
            set(get(AX(2),'Ylabel'),'String','d(spot #))') 
            saveas(1,[im_folder,sub_folder,imname,foci_add,th_add,fig_tail]);
        
        figure(2)
            clf
            subplot(1,2,1)
            bar(xout0,n0)
            title(im_folder)
            xlabel('Intensity (A.U.)')
            ylabel('#')
            set(gca,'YScale','log')
            %%%%%%%%%%%%%%%%%
            subplot(1,2,2)
            bar(xout,n)
            title(imname,'Interpreter','none')
            xlabel('Intensity (A.U.)')
            ylabel('#')
            set(gca,'YScale','log')
            saveas(2,[im_folder,sub_folder,imname,fig_tail]);
        
        data_name = [imname,data_tail];
        try
            delete([im_folder,sub_folder,data_name])
        catch
        end
        xlswrite([im_folder,sub_folder,data_name],max_spot);
    end
end