function [lld,varargout] = FSP_NSX_codegrade_total_P0(k,Nstate,ndata,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to calculate the likelihood of total mRNAs from a 2-state-elongation model with co-degradation: %%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% k (input): [transition rates, degradation rate, P0 extra, t2, t3]: the transition rates can be 2 forms:
%%              1. transition matrix: kk_ij = k_ij (i~=j), kk_ii = kTX_i;
%%              2. free-ordered elements: the indices in the transition matrix is defined by Lk.
%% Nstate (input): number of states.
%% ndata (input): experimental data (nascent mRNA#)
%% varargin (input): 1. tel2: post-signal elongation time;
%%                   2. tel3: pre-signal elongation time;
%%                   3. Lk: indices of kk in the transition matrix.
%%                   4. pTr: FISH probe binomial binding matrix.
%%                   5. Obin: binning for llp calculation.
%%                   6. TXmax: Max of mRNA # for llp calculation.
%%                   7. pc: distribution of gene copy number.
%%                   8. tr: time after rifampicin application. 
%% lld (output): log-likelihood of the experimental data.
%% varargout (output): 1. [nRNA1',PP0']: nascent mRNA # (bin to 1) and it probability distribution;
%%                     2. [nRNA',reshape(PP,Nstate,Nx)']: nascent mRNA # (finer binning) and the probability distribution considering different gene states.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k(k > 5000) = 5000;

if isempty(varargin) || length(varargin) <= 2 || isempty(varargin{3})
    Lk = [2:(Nstate+1):Nstate^2,(Nstate+1):(Nstate+1):Nstate^2,1:(Nstate+1):Nstate^2];
else
    Lk = varargin{3};
end
if length(varargin) <= 3 || isempty(varargin{4})
    pTr = [];
    Nmax0 = 100;
else
    pTr = varargin{4};
    Nmax0 = size(pTr,1)-1;
end

if length(varargin) <= 4 || isempty(varargin{5})
    Obin = 10;
else
    Obin = varargin{5};
end

if length(varargin) <= 5 || isempty(varargin{6})
    TXmin = 0;
    TXmax = 100;
elseif length(varargin{6}) == 1
    TXmin = 0;
    TXmax = varargin{6};
else
    TXmin = varargin{6}(1);
    TXmax = varargin{6}(2);
end

Nmax0 = max(Nmax0,TXmax);

if length(varargin) <= 6 || isempty(varargin{7})
    pc = 1;
else
    pc = varargin{7};
end

if length(varargin) <= 7 || isempty(varargin{8})
    tr = 0;
else
    tr = varargin{8};
end

Lt = ismember(Lk,[1:(Nstate+1):Nstate^2]);
k00 = k(1:numel(Lk));
ks = k00(~Lt);
kr = zeros(1,Nstate);
kr(ceil(Lk(Lt)/Nstate)) = k00(Lt);   %%% kr1,kr2,kr3,...

NL = numel(Lk);
kd = k(NL+1);

if numel(k) > NL+1
    P0e = k(NL+2);
else
    P0e = 0;
end

if P0e > 1 || P0e < 0
    error(['Out of range: P0_extra = ',num2str(P0e)])
end

if (isempty(varargin) || isempty(varargin{1})) && numel(k) > NL+2
    t2 = k(NL+3);
elseif ~isempty(varargin) && ~isempty(varargin{1})
    t2 = varargin{1};
else
    t2 = 0;
end

if (numel(varargin) <= 1 || isempty(varargin{2})) && numel(k) > NL+3
    t3 = k(NL+4);
elseif numel(varargin) > 1 && ~isempty(varargin{2})
    t3 = varargin{2};
else
    t3 = 0;
end

Km = zeros(Nstate);
Km(Lk(~Lt)) = ks;
Km = Km-diag(sum(Km));
Km0 = Km;
Km0(1,:) = 1;
Pz = zeros(Nstate,1);
Pz(1) = 1;
P0 = Km0\Pz;   %%% equilibrium state distribution

Nmax = Nmax0;   %%% maximal number of complete RNA molecules
Nbin = 10;   %%% Number of RNA bin for 1 mRNA;
if max(k) <= 200
    Nbint = 1000;   %%% Number of time bin
else
    Nbint = ceil(max(k)/1000)*10000;
end
Nbint2 = floor(Nbint*t2);   %%% Number of time bin in t2 (post_label)
Nbin2 = floor(Nbin*t2);   %%% Number of RNA bin in t2 (post_label)
% dt2 = t2-Nbint2/Nbint;
Nbint3 = floor(Nbint*t3);   %%% Number of time bin in t3 (pre_label)
% dt3 = t3-Nbint3/Nbint;
Nr = Nbint/Nbin;
dt = 1/Nbint;
dRNA = 1/Nbin;
nRNA = 0:dRNA:Nmax;   %%% number of complete RNA molecules
Nx = length(nRNA);


KQ = Km*dt+diag(-kr*dt+1);
Q0 = blkdiag0(KQ,Nx,Nstate);
% Q0 = kron(sparse(1:Nx,1:Nx,1),KQ);

% temp = cell(1,Nx);
% [temp{:}] = deal(sparse(KQ));
% Q0 = blkdiag(temp{:});

LRNA = [0:dRNA:1]';
deg = [exp(-kd*LRNA);0];

LRNA2 = [0:dRNA:t2]';
deg2 = exp(-kd*LRNA2);
d0 = exp(kd*dRNA);

ddeg0 = -diff(deg);
ddeg = ddeg0*deg2(end);
ddeg(1) = ddeg(1)+1-deg2(end);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initial condition setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % KQ1 = Km+diag(-kr-kd);
% % Q01 = blkdiag0(KQ1,Nx,Nstate);
% % ddeg21 = kron(ddeg,diag(kr));
% % dQ01 = blkdiag2(ddeg21,Nx,Nstate);
% % ddeg22 = kron(ddeg,diag(kd*ones(size(kr))));
% % dQ02 = transpose(blkdiag2(ddeg22,Nx,Nstate));
% % QQ1 = Q01+dQ01+dQ02;
% % 
% % QQ1(1,:) = 1;
% % Pz = sparse(1,1,1,Nx*Nstate,1);
% % PP = full(QQ1\Pz);   %%% equilibrium state distribution

Nk4 = max(4,ceil((tr-1-t2-t3)*kd)+1);

t4 = Nk4/kd;
Nbint4 = floor(Nbint*t4);   %%% Number of time bin in t4 (mature time)
Nbin4 = floor(Nbin*t4);   %%% Number of RNA bin in t4 (mature time)
LRNA4 = [0:dRNA:t4]';
deg4 = exp(-kd*LRNA4);

PP = zeros(Nx*Nstate,1);
PP(1:Nstate) = P0;   %%% initial state

if tr > 0
    Nbintr = floor(Nbint*(1+t2+t3+t4-tr));   %%% Number of time bin for actual propagation (considering rifampicin appication)
else
    Nbintr = Nbint4+Nbint2+Nbint+Nbint3;
end

iel = Nbin4;
iel0 = iel;
ddeg2 = kron(ddeg,diag(kr*dt));
dQ0 = blkdiag2(ddeg2,Nx,Nstate);
M0 = sparse(1:Nx*Nstate,1:Nx*Nstate,repmat(kr*dt,1,Nx));
QQ = Q0+dQ0*deg4(iel+1)+M0*(1-deg4(iel+1));
for ii = Nbint4:-1:max(1,(Nbint4-Nbintr+1))
    iel = round((ii-1)/Nr);
    if iel < iel0
        QQ = Q0+dQ0*deg4(iel+1)+M0*(1-deg4(iel+1));
        iel0 = iel;
    end
    PP = QQ*PP;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Propagation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t2 > 0
    iel = Nbin2;
    iel0 = iel;
    ddeg2 = kron(ddeg,diag(kr*dt));
    dQ0 = blkdiag2(ddeg2,Nx,Nstate);
    M0 = sparse(1:Nx*Nstate,1:Nx*Nstate,repmat(kr*dt*(1-d0),1,Nx));
    QQ = Q0+dQ0;
    for ii = Nbint2:-1:max(1,(Nbint4+Nbint2-Nbintr+1))
        iel = round((ii-1)/Nr);
        if iel ~= iel0
            dQ0 = dQ0*d0+M0;
            QQ = Q0+dQ0;
            iel0 = iel;
        end
        PP = QQ*PP;
    end
end

iel = Nbin;
iel0 = iel;
ddeg2 = kron(ddeg0,diag(kr*dt));
dQ0 = blkdiag2(ddeg2,Nx,Nstate);
QQ = Q0+dQ0;
for ii = Nbint:-1:max(1,(Nbint4+Nbint2+Nbint-Nbintr+1))
    iel = round((ii-1)/Nr);
    if iel ~= iel0
        temp = kr'*dt*deg(iel+2);
        temp1 = temp(:,ones(1,Nx-iel));
        temp1 = temp1(:)';
        QQ([(iel0*Nstate+1):Nx*Nstate]+[0:((Nx-iel0)*Nstate-1)]*Nx*Nstate) = 0;
        QQ([(iel*Nstate+1):Nx*Nstate]+[0:((Nx-iel)*Nstate-1)]*Nx*Nstate) = QQ([(iel*Nstate+1):Nx*Nstate]+[0:((Nx-iel)*Nstate-1)]*Nx*Nstate)+temp1;
        iel0 = iel;
    end
    PP = QQ*PP;
end

if Nbint3 > 0
    for ii = Nbint3:-1:max(1,(Nbint4+Nbint2+Nbint+Nbint3-Nbintr+1))
        PP = QQ*PP;
    end
end

% PP(PP<0) = 0;
% PP = PP/sum(PP);
tfalse = any(PP < 0);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% log likelihood calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PP = PP*(1-P0e);
PP(1) = PP(1)+P0e;

PP0 = sum(reshape(PP,Nstate,Nx));

%%% consider multiple gene copies:
if length(pc) > 1
    ptemp = PP0;
    pall = PP0*pc(1);
    for ii = 2:length(pc)
        ptemp = conv(ptemp,PP0);
        pall = ptemp*pc(ii)+[pall,zeros(1,size(PP0,2)-1)];
    end
    PP0 = pall(1:size(PP0,2));
end

PP0 = conv(PP0,ones(1,Obin));
% PP0 = (PP0(round(Nbin/2):Nbin:end)+PP0(round(Nbin/2)+1:Nbin:end))/2;
PP0 = (PP0(round(Obin/2):Obin:end)+PP0(round(Obin/2)+1:Obin:end))/2;
if ~isempty(pTr)
    PP0 = PP0*pTr;
end

nRNA1 = nRNA(1:Obin:end);
% Nbin = 1;

logPP = log(PP0);
% NnRNA = hist(ndata,nRNA);
if tfalse
    lld = inf;
else
%     lld = -NnRNA*logPP';
% %     ndata(ndata > Nmax) = Nmax;
%     ndata(ndata > TXmax) = TXmax;
    ndata(ndata < 0) = 0;
    ndata = ndata((ndata <= TXmax) & (ndata >= TXmin));
    lld = -sum(logPP(round(ndata*Nbin/Obin)+1));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargout = {[nRNA1',PP0'],[nRNA',reshape(PP,Nstate,Nx)']};



function Q0 = blkdiag0(KQ,Nx,Nstate)
% [Ix0,Iy0] = ind2sub([Nstate,Nstate],[1:Nstate^2]');
Ix0 = repmat([1:Nstate]',1,Nstate);
Iy0 = repmat([1:Nstate],Nstate,1);
% % [Iy0,Ix0] = meshgrid(1:Nstate);
Ix = repmat(Ix0(:),1,Nx)+repmat([0:Nstate:Nstate*(Nx-1)],Nstate*Nstate,1);
Iy = repmat(Iy0(:),1,Nx)+repmat([0:Nstate:Nstate*(Nx-1)],Nstate*Nstate,1);
K0 = repmat(KQ(:),1,Nx);
Q0 = sparse(Ix(:),Iy(:),K0(:),Nx*Nstate,Nx*Nstate);



function dQ0 = blkdiag2(ddeg2,Nx,Nstate)
Ix0 = bsxfun(@plus,[1:size(ddeg2,1)]',kron([0:(Nx-1)]*Nstate,ones(1,Nstate)));
Iy0 = repmat(1:(Nx*Nstate),size(ddeg2,1),1);
K0 = repmat(ddeg2,1,Nx);

Ix = Ix0(:);
Iy = Iy0(:);
K0 = K0(:);
Itrue = Ix <= Nx*Nstate;

dQ0 = sparse(Ix(Itrue),Iy(Itrue),K0(Itrue),Nx*Nstate,Nx*Nstate);

