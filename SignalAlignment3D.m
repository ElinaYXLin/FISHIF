function [R_m,T_m]=SignalAlignment3D(Benchmark_data,Calibration_data)
%% Program Intro
% a program to Match data.
% SignalAlignment3D(Benchmark_data,Calibration_data)
%% Test & correct data
if size(Benchmark_data,2)>3
Benchmark_data=Benchmark_data';
end
if size(Calibration_data,2)>3
Calibration_data=Calibration_data';
end
%% centre
% centre_b=mean(Benchmark_data,2);
% centre_c=mean(Calibration_data,2);
% B_centre=Benchmark_data-centre_b;
% C_centre=Calibration_data-centre_c;
%% Coherent point drift
fixed=Benchmark_data;
moving=Calibration_data;
%option
opt.method='rigid';
opt.scale=0;
opt.tol=1e-15;
opt.max_it=1e2;
Transform=cpd_register(fixed,moving,opt);
Transform.X=fixed;
R_m=Transform.R;
T_m=Transform.t;
%% visible
% Initial point-sets
figure,cpd_plot_iter(fixed, moving); title('Before');
% Registered point-sets
figure
cpd_plot_iter(fixed, Transform.Y);
title('After');
legend('Fixed point cloud','Moving point cloud')
end