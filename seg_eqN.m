function [nucleus_bin0,bin_min0,bin_max0] = seg_eqN(nucleus_distance,Nseg,average_radius)

nucleus_distance = sort(nucleus_distance);
Nall = length(nucleus_distance);
Nr = round(Nall*average_radius);

Nc = (Nr+1):((Nall-2*Nr-1)/(Nseg-1)):(Nall-Nr);
N_min = round(Nc-Nr);
N_max = round(Nc+Nr);
N_min2 = max(1,round(Nc-Nr)-1);
N_max2 = min(Nall,round(Nc+Nr)+1);

bin_min0 = (nucleus_distance(N_min)+nucleus_distance(N_min2))/2;
bin_max0 = (nucleus_distance(N_max)+nucleus_distance(N_max2))/2;
nucleus_bin0 = zeros(1,Nseg);
for ii = 1:Nseg
    nucleus_bin0(ii) = mean(nucleus_distance(N_min(ii):N_max(ii)));
end
