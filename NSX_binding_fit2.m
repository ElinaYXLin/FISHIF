function [kk,varargout] = NSX_binding_fit2(xdata,ydata,kk0,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to fit data to (N+1)-state TF binding model %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% xdata (input): TF concentration
%% ydata (input): Occupancies of different binding states (conentration,binding states) (binding states from 0 to N)
% % %% kk0 (input): initial values of fitting parameters ((c0,h0,ki) for binding states 1 to N, kn0)
%% kk0 (input): initial values of fitting parameters ((c0,h0,ki) for binding states 1 to N, and k20 (optional),ki0)
%% varargin (input): 1. number of trials (default: 10); 2. dI
%% kk (output): Results of fitting parameters
%% varargout(output): {kk_all,fval_all,exitflag_all}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % lb = zeros(size(kk0)); %lb(3:3:end-2) = kk0(3:3:end-2);% lb(end-1:end) = 0;
% % % lb = zeros(size(kk0)); lb(3:3:end-2) = kk0(3:3:end-2);% lb(end-1:end) = 0;
% % % ub = [kk0(1:end-2)*3,20,20]; ub(2:3:end-2) = kk0(2:3:end-2)*3; ub(3:3:end-2) = kk0(3:3:end-2);
% % ub = [kk0(1:end-2)*3,20]; ub(2:3:end-1) = 7; ub(3:3:end-1) = 20;
% % % ub = [kk0(1:end-2)*3,20,20]; ub(2:3:end-2) = 7; ub(3:3:end-2) = kk0(3:3:end-2);
epsilon = 1e-6;

if isempty(varargin)
    Nfit = 10;
else
    Nfit = varargin{1};
end

if length(varargin) < 2 || isempty(varargin{2})
    dI = mod(length(kk0),3);
else
    dI = varargin{2};
end

c0 = kk0(1:3:end-dI); c0l = zeros(size(c0))+epsilon; c0u = c0*3;
h0 = kk0(2:3:end-dI); h0l = h0; h0u = h0;
% h0 = kk0(2:3:end-dI); h0l = h0*0; h0u = h0*5;
ki = kk0(3:3:end-dI); kil = zeros(size(ki))+epsilon; kiu = 10*ones(size(ki));

% % if dI == 2
% %     k20 = kk0(end-1); k20l = 0+epsilon; k20u = 20;
% %     ki0 = kk0(end); ki0l = 0+epsilon; ki0u = 20;
% % elseif dI == 1
% %     k20 = []; k20l = []; k20u = [];
% %     ki0 = kk0(end); ki0l = 0+epsilon; ki0u = 20;
% % end

ki0 = kk0(end); ki0l = 0+epsilon; ki0u = 20;
% ks0 = kk0(end-dI+1:end-1); ks0l = zeros(size(ks0)); ks0u = zeros(size(ks0));
ks0 = kk0(end-dI+1:end-1); ks0l = zeros(size(ks0)); ks0u = 20*ones(size(ks0));

lb = zeros(size(kk0)); ub = zeros(size(kk0)); 
lb(1:3:end-dI) = c0l; ub(1:3:end-dI) = c0u; 
lb(2:3:end-dI) = h0l; ub(2:3:end-dI) = h0u; 
% lb(2:3:end-dI) = zeros(size(h0l)); ub(2:3:end-dI) = 2*h0u; 
lb(3:3:end-dI) = kil; ub(3:3:end-dI) = kiu; 
% % if dI == 2
% %     lb(end-1) = k20l; ub(end-1) = k20u;
% % end
lb(end-dI+1:end) = [ks0l,ki0l]; ub(end-dI+1:end) = [ks0u,ki0u];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% MLE Fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk_all = zeros(Nfit,length(kk0));
fval_all = zeros(Nfit,1);
exitflag_all = zeros(Nfit,1);

% options = SIMPSASET('COOL_RATE',0.5,'TEMP_END',1);
options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',15000000,'FunctionTolerance',1e-10);
% % [kk,FVAL,EXITFLAG] = simulannealbnd(@(x) -sum(sum(log(NSX_binding(x,xdata)).*ydata)),kk0,lb,ub,options);

tic
pool_name = parpool(min(30,Nfit));
parfor ii = 1:Nfit
    [kk_all(ii,:),fval_all(ii),exitflag_all(ii)] = simulannealbnd(@(x) log_NSX_bind(x,xdata,ydata,dI),kk0,lb,ub,options);
%     disp(['N = ',num2str(ii),'/',num2str(Nfit),', t = ',num2str(toc)])
end
delete(pool_name)
toc

[~,Imin] = min(fval_all);
kk = kk_all(Imin,:);

varargout = {kk_all,fval_all,exitflag_all};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc


function lld = log_NSX_bind(kk,cx,Nx,dI)
% lld = -sum(sum(log(NSX_binding(kk,cx)).*Nx));
lld = -sum(sum(log(NSX_binding2(kk,cx,dI)).*Nx));



