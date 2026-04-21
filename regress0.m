function [kfit,dkfit] = regress0(y0,x0,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to calculate the optimal range of linear regression of data %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kfit (output): regression parameters.
%% dkfit (output): upper and lower bounds of regression parameters within certain range of the squared error.
%% y0 (input): input y for regression.
%% x0 (input): input x for regression.
%% varargin (input): {epsilon} relative range of the squared error with respect to the minimum value.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin) || isempty(varargin{1})
    epsilon = 0.05;
else
    epsilon = varargin{1};
end

k0 = (x0'*x0)\(x0'*y0);

Smin = (y0-x0*k0)'*(y0-x0*k0);
[E0,L0] = eig(x0'*x0);
dk2 = sqrt(epsilon*Smin./diag(L0));
dk = E0*diag(dk2);
dk0 = max(abs(dk),[],2);

kfit = k0;
dkfit = [k0-dk0, k0+dk0];



