clear all
close all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_seg = 1000;   %%% Number of segment in a nucleus
N_nu = 1000;   %%% Number of simulated nuclei
M_max = 10000;   %%% Maximal number of particles of a nucleus
photon = 100;   %%% Mean number of photon emitted per protein
M_max2 = 800;   %%% Maximal number of particle2 of a nucleus
photon2 = 180;   %%% Mean number of photon emitted per protein2
r_decay = 0.2;   %%% Decay length of protein profile
N_sec = 5;   %%% Number of secondary antibody binding sites per primary antibody
p_sec = 0.8;   %%% probability of secondary antibody binding
p0 = 0.1;   %%% base probability of secondary antibody binding
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_num = zeros(N_nu,1);
var_num = zeros(N_nu,1);
EL = zeros(N_nu,1);
mean_num0 = zeros(N_nu,1);

for I_nu = 1:N_nu
    EL(I_nu) = rand(1);
    M_protein = ceil(M_max*(exp(-EL(I_nu)/r_decay)-exp(-1./r_decay))/(1-exp(-1./r_decay)));
    position_protein = ceil(N_seg*rand(M_protein,1));
    N_position =  poissrnd(photon*binornd(N_sec*hist(position_protein,[1:N_seg]),p_sec*((1-p0)*exp((EL(I_nu)-1)/r_decay)+p0)));
%     N_position =  poissrnd(photon*binornd(N_sec*hist(position_protein,[1:N_seg]),p_sec));

    M_protein2 = ceil(M_max2*rand(1));
    position_protein2 = ceil(N_seg*rand(M_protein2,1));
%     N_position2 =  poissrnd(photon2*binornd(N_sec*hist(position_protein2,[1:N_seg]),p_sec*((1-p0)*exp((1-EL(I_nu))/r_decay)+p0)));
    N_position2 =  poissrnd(photon2*binornd(N_sec*hist(position_protein2,[1:N_seg]),p_sec));
    
    mean_num(I_nu) = mean(N_position);
    var_num(I_nu) = var(N_position);
    mean_num0(I_nu) = mean(N_position);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(mean_num,var_num,'.')
xlabel('Mean intensity (A.U.)')
ylabel('Intensity variance (A.U.^2)')
title('Simulated nuclear intensity plot')

% p = polyfit(mean_num,var_num,1);
% figure(2)
% plot(mean_num,(var_num-p(1)*mean_num-p(2)).^2,'.')

p_post = polyfit(mean_num(EL > 0.8),var_num(EL > 0.8),1);
p_ant = polyfit(mean_num(EL < 0.5),var_num(EL < 0.5),1);
% mean_num_re = (mean_num*p_post(1)-var_num+p_post(2))/(p_post(1)-p_ant(1));
mean_num_re = (mean_num*p_post(1)-var_num)/(p_post(1)-p_ant(1));
p_re = polyfit(mean_num0,mean_num_re,1);
x_range = [0:max(mean_num0)/1000:max(mean_num0)];

figure(2)
plot(mean_num0,mean_num_re,'b.',x_range,p_re(1)*x_range+p_re(2),'r-')
legend('raw data',['Linear fit: y = ',num2str(p_re(1)),'*x + ',num2str(p_re(2))])
xlabel('Mean intensity (A.U.)')
ylabel('Reconstructed intensity (A.U.)')
title('Intensity reconstruction')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

