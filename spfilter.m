function [max_spot,varargout] = spfilter(peak_parameter,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% A function to filter out low intensity/large eccentricity noise: %%%
%% max_spot: spot characteristics (Intensity,sigma_x,sigma_y,sigma_z)
%% varargout (output): {I_real} indices for real spots
%% peak_parameter: fitting parameter values for each local maxima with the 
%%                 form: (I0,sigma_x,sigma_y,sigma_z,A,x0,y0,z0,theta,residual)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmin = 0.4;
rzmin = 0.25;

true_spot = true(size(peak_parameter,1),1);
true_spot = true_spot & (peak_parameter(:,1) <= 1e5) & (peak_parameter(:,1) > 0);
true_spot = true_spot & (peak_parameter(:,2) > 0) & (peak_parameter(:,3) > 0) & (peak_parameter(:,4) > 0);

if length(varargin) >= 4 && ~isempty(varargin{4})
    true_spot = true_spot & (prod(peak_parameter(:,1:3),2)*2*pi >= varargin{4});
end

if length(varargin) >= 3 && ~isempty(varargin{3})
    true_spot = true_spot & (sqrt(1-(peak_parameter(:,2)./peak_parameter(:,3)).^2) <= varargin{3});
else
    true_spot = true_spot & (sqrt(1-(peak_parameter(:,2)./peak_parameter(:,3)).^2) <= 0.9);
end

if isempty(varargin)
    true_spot = true_spot & (peak_parameter(:,2) >= rmin) & (peak_parameter(:,2) <= 3);
    true_spot = true_spot & (peak_parameter(:,3) >= rmin) & (peak_parameter(:,3) <= 3);
elseif length(varargin{1}) == 1 
    true_spot = true_spot & (peak_parameter(:,2) >= rmin) & (peak_parameter(:,2) <= varargin{1});
    true_spot = true_spot & (peak_parameter(:,3) >= rmin) & (peak_parameter(:,3) <= varargin{2});
else
    rxy = sqrt(peak_parameter(:,2).*peak_parameter(:,3));
    true_spot = true_spot & (rxy >= varargin{1}(1)) & (rxy <= varargin{1}(2));
end

if length(varargin) >= 5 && ~isempty(varargin{5})
    true_spot = true_spot & (peak_parameter(:,4) >= varargin{5}(1)) & (peak_parameter(:,4) <= varargin{5}(2));
else
    true_spot = true_spot & (peak_parameter(:,4) >= rzmin) & (peak_parameter(:,4) <= 2);
end
%true_spot = true_spot & (peak_parameter(:,4) > 0.25) & (peak_parameter(:,3) <2);


%max_spot(:,1) = peak_parameter(true_spot,1);
%max_spot(:,1) = peak_parameter(:,1).*peak_parameter(:,2).*peak_parameter(:,3)*sqrt(2*pi);
%max_spot(:,2:10) = peak_parameter(true_spot,2:10);
max_spot = peak_parameter(true_spot,:);

varargout = {true_spot};
