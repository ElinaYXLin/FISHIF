function [pp,varargout] = NSX_binding2(kk,cx,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to calculate the occupancy at certain TF concentration for (N+1)-state TF binding model %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% kk (input): model parameters ((c0,h0,ki) for binding states 1 to N, kn0, and a0 for the kn0/k01)
%% kk (input): model parameters ((c0,h0,ki) for binding states 1 to N, and k20 (optional),ki0)
%% cx (input): TF concentration
%% varargin (input): {dI} 
%% pp (output): Binding occupancies for different TF binding states (cx,pp)
%% varargout (output): Binding occupancies for different TF binding states (pp0, including hidden state)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin) || isempty(varargin{1})
    dI = 2;
else
    dI = varargin{1};
end

%% 
c00 = kk(1:3:end-dI);
h00 = kk(2:3:end-dI);
ki0 = kk(3:3:end-dI);
kn0 = kk(end)*ones(size(ki0)); 
if dI >= 2
    kn0(1:dI-1) = kk(end-dI+1:end-1);
end
k01 = 1;

cc = repmat(cx,1,numel(c00));
c0 = repmat(c00,numel(cx),1);
h0 = repmat(h00,numel(cx),1);
ki00 = repmat(ki0,numel(cx),1);
kk2 = ki00.*(cc./c0).^h0;

Nk0 = numel(c00)+2;
pp0 = zeros((size(cc)+[0,2]));
A0 = zeros(Nk0,1); A0(1) = 1;
for ii = 1:numel(cx)
    kk_all = zeros(Nk0);
    kk_all(2:(Nk0+1):end) = [kk2(ii,:),0];
    kk_all((Nk0+1):(Nk0+1):end-2) = ki0;
    kk_all(Nk0,2:end-1) = kn0;
    kk_all(1,end) = k01;
    
    kk_all(1:(Nk0+1):end) = -sum(kk_all);
    kk_all(1,:) = 1;
    
    pp0(ii,:) = kk_all\A0;
end

pp = pp0(:,1:end-1);
pp(:,1) = pp0(:,1)+pp0(:,end);

varargout{1} = pp0;

