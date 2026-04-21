function color2cmp(lf_name)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to compare two color foci quantification results: %%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
close all

sub_folder = 'Histogram/';
sub_folder2 = 'Histogram2/';
sub_folderA = 'Histogram_A/';
sub_folder2A = 'Histogram2_A/';

in_folder = 'stacks/';
input_name = 'matchlist.xls';
data_tail = '_raw.xls';
output_add = '_spot_fit';
mat_tail = '.mat';
fig_tail = '.fig';
dmax = 2;   %%% maximal xy mismatch

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
        data_name = [sub_list{list_J,3}(1:end-1),data_tail];
        raw_data = xlsread([im_folder,sub_folder,data_name]);
        Inten_spot1 = raw_data(:,1).*raw_data(:,2).*raw_data(:,3)*2*pi;
        xyz1 = raw_data(:,6:8);
        
        data_name = [sub_list{list_J,3}(1:end-1),data_tail];
        raw_data = xlsread([im_folder,sub_folder2,data_name]);
        Inten_spot2 = raw_data(:,1).*raw_data(:,2).*raw_data(:,3)*2*pi;
        xyz2 = raw_data(:,6:8);
        
        load([im_folder,sub_folderA,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'b')
        Inten_spot1 = Inten_spot1/b;
        load([im_folder,sub_folder2A,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'b')
        Inten_spot2 = Inten_spot2/b;
        
        %%% Spot matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dxy = pdist2(xyz1(:,1:2),xyz2(:,1:2));   %%% xy distance of RNA spots in two channels
        I_match = zeros(0);
        for I1 = 1:size(xyz1,1)
            [dmin,Imin] = min(dxy(I1,:));
            if dmin <= dmax && Inten_spot1(I1) <= 40 && Inten_spot2(Imin) <= 40 
                I_match = [I_match;[I1,Imin]];
            end
        end
        
        Inten_spot1 = Inten_spot1(I_match(:,1));
        Inten_spot2 = Inten_spot2(I_match(:,2));
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [sum(Inten_spot1),sum(Inten_spot2)]
        figure
            plot(Inten_spot1,Inten_spot2,'k.')
            hold on
            p = polyfit(Inten_spot1,Inten_spot2,1);
%             p1 = cov(Inten_spot1,Inten_spot2)/var(Inten_spot1);
%             p = [p1(1,2),0];
            plot(Inten_spot1,Inten_spot1*p(1)+p(2),'r')
            legend('Foci',['Fit: k = ',num2str(p(1))])
            xlabel('Channel A (#)')
            ylabel('Channel B (#)')
% %         saveas(gcf,[im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,fig_tail])
        
    end
end