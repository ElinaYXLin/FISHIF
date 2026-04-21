function [kk,varargout] = NSX_binding_fit2b(xdata,ydata,kk0,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to fit data to (N+1)-state TF binding model %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% xdata (input): TF concentration for individual nuclei
%% ydata (input): Occupancies of different binding states (conentration,binding states) (binding states from 0 to N)
% % %% kk0 (input): initial values of fitting parameters ((c0,h0,ki) for binding states 1 to N, kn0)
%% kk0 (input): initial values of fitting parameters ((c0,h0,ki) for binding states 1 to N, and k20,ki0)
%% varargin (input): number of trials (default: 10)
%% kk (output): Results of fitting parameters
%% varargout(output): {kk_all,fval_all,exitflag_all}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % lb = zeros(size(kk0)); %lb(3:3:end-2) = kk0(3:3:end-2);% lb(end-1:end) = 0;
% % % lb = zeros(size(kk0)); lb(3:3:end-2) = kk0(3:3:end-2);% lb(end-1:end) = 0;
% % % ub = [kk0(1:end-2)*3,20,20]; ub(2:3:end-2) = kk0(2:3:end-2)*3; ub(3:3:end-2) = kk0(3:3:end-2);
% % ub = [kk0(1:end-2)*3,20]; ub(2:3:end-1) = 7; ub(3:3:end-1) = 20;
% % % ub = [kk0(1:end-2)*3,20,20]; ub(2:3:end-2) = 7; ub(3:3:end-2) = kk0(3:3:end-2);

c0 = kk0(1:3:end-2); c0l = zeros(size(c0)); c0u = c0*3;
h0 = kk0(2:3:end-2); h0l = h0; h0u = h0;
ki = kk0(3:3:end-2); kil = zeros(size(ki)); kiu = 20*ones(size(ki));
k20 = kk0(end-1); k20l = 0; k20u = 20;
% % k20 = kk0(end-1); k20l = 0; k20u = 0;
ki0 = kk0(end); ki0l = 0; ki0u = 20;

lb = zeros(size(kk0)); ub = zeros(size(kk0)); 
lb(1:3:end-2) = c0l; ub(1:3:end-2) = c0u; 
% lb(2:3:end-1) = h0l; ub(2:3:end-2) = h0u; 
lb(2:3:end-2) = zeros(size(h0l)); ub(2:3:end-2) = 2*h0u; 
lb(3:3:end-2) = kil; ub(3:3:end-2) = kiu; 
lb(end-1) = k20l; ub(end-1) = k20u;
lb(end) = ki0l; ub(end) = ki0u;

if isempty(varargin)
    Nfit = 10;
else
    Nfit = varargin{1};
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% MLE Fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk_all = zeros(Nfit,length(kk0));
fval_all = zeros(Nfit,1);
exitflag_all = zeros(Nfit,1);

% options = SIMPSASET('COOL_RATE',0.5,'TEMP_END',1);
options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',15000000,'FunctionTolerance',1e-10);
% % [kk,FVAL,EXITFLAG] = simulannealbnd(@(x) -sum(sum(log(NSX_binding(x,xdata)).*ydata)),kk0,lb,ub,options);

tic
for ii = 1:Nfit
    [kk_all(ii,:),fval_all(ii),exitflag_all(ii)] = simulannealbnd(@(x) log_NSX_bind(x,xdata,ydata),kk0,lb,ub,options);
    disp(['N = ',num2str(ii),'/',num2str(Nfit),', t = ',num2str(toc)])
end

[~,Imin] = min(fval_all);
kk = kk_all(Imin,:);

varargout = {kk_all,fval_all,exitflag_all};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc


function lld = log_NSX_bind(kk,cx,Nx)
% % % lld = -sum(sum(log(NSX_binding(kk,cx)).*Nx));
% % lld = -sum(sum(log(NSX_binding2(kk,cx)).*Nx));
p00 = zeros(size(Nx));
for ii = 1:length(cx)
    p00(ii,:) = mean(NSX_binding2(kk,cx{ii}));
end
lld = -sum(sum(log(p00).*Nx));


