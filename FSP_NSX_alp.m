function [lld,varargout] = FSP_NSX_alp(k,Nstate,ndata,t2,t3,varargin)
%FSP_NSX(x,Ns,f_ob,t2,[],Lk)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to calculate the likelihood of 2-state-elongation model: %%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% k (input): [transition rates, degradation rate]: the transition rates can be 2 forms:
%%              1. transition matrix: kk_ij = k_ij (i~=j), kk_ii = kTX_i;
%%              2. free-ordered elements: the indices in the transition matrix is defined by Lk.
%% Nstate (input): number of states.
%% ndata (input): experimental data (nascent mRNA#)
%% varargin (input): 1. tel2: post-signal elongation time;
%%                   2. tel3: pre-signal elongation time;
%%                   3. Lk: indices of kk in the transition matrix.
%%                   4. pTr or [pTr, TXmax]: PTr: FISH probe binomial binding matrix. TXmax:range (or max) of mRNA # for llp calculation.
%%                   5. Nth: nascent mRNA detection threshold.
%% lld (output): log-likelihood of the experimental data.
%% varargout (output): 1. [nRNA1',PP0']: nascent mRNA # (bin to 1) and it probability distribution;
%%                     2. [nRNA',reshape(PP,Nstate,Nx)']: nascent mRNA # (finer binning) and the probability distribution considering different gene states.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k(k > 5000) = 5000;
alp=k(end);k=k(1:end-1);
if isempty(varargin) || length(varargin) <= 2 || isempty(varargin{3})
    Lk = [2:(Nstate+1):Nstate^2,(Nstate+1):(Nstate+1):Nstate^2,1:(Nstate+1):Nstate^2];
else
    Lk = varargin{3};
end
if isempty(varargin) || length(varargin) <= 3 || isempty(varargin{4})
    pTr = [];
    Nmax0 = 100;
    TXmax = 100;
else
    if length(varargin{4}) == 1
        pTr = varargin{4};
        Nmax0 = size(pTr,1)-1;
        TXmax = 60;
    elseif ~isempty(varargin{4}{1})
        pTr = varargin{4}{1};
        Nmax0 = size(pTr,1)-1;
        TXmax = varargin{4}{2};
    else
        pTr = [];
        Nmax0 = varargin{4}{2};
        TXmax = varargin{4}{2};
    end
end

if isempty(varargin) || length(varargin) <= 4
    Nth = 4;
else
    Nth = varargin{5};
end
%Nstate=2
Lt = ismember(Lk,[1:(Nstate+1):Nstate^2]);
k00 = k(1:numel(Lk));
ks = k00(~Lt);%k01?k10
kr = zeros(1,Nstate);
kr(ceil(Lk(Lt)/Nstate)) = k00(Lt);   %%% kr1,kr2,kr3,...kini1,kini2

NL = numel(Lk);
%if (isempty(varargin) || isempty(varargin{1})) && numel(k) > NL
%    t2 = k(NL+1);
%elseif ~isempty(varargin) && ~isempty(varargin{1})
%    t2 = varargin{1};
%else
%    t2 = 0;
%end
%t2=0.4
%if (numel(varargin) <= 1 || isempty(varargin{2})) && numel(k) > NL+1
%    t3 = k(NL+2);
%elseif numel(varargin) > 1 && ~isempty(varargin{2})
%    t3 = varargin{2};
%else
%    t3 = 0;
%end

Km = zeros(Nstate);
Km(Lk(~Lt)) = ks;
Km = Km-diag(sum(Km));
Km0 = Km;
Km0(1,:) = 1;
Pz = zeros(Nstate,1);
Pz(1) = 1;
P0 = Km0\Pz;   %%% equilibrium state distribution%%%??
% P0(1)=P0(1)+randn*P0(1);
% P0(2)=1-P0(1);
% P0(P0<0)=0;
% P0=p_ini;
Nmax = Nmax0;   %%% maximal number of complete RNA molecules
Nbin = 10;   %%% Number of RNA bin for 1 mRNA;%%%%%%???1???MRNA????
if max(k) <= 200
    Nbint = 1000;   %%% Number of time bin
else
    Nbint = ceil(max(k)/1000)*10000;
end
Nbint = 2000;
Nbint2 = floor(Nbint*t2);   %%% Number of time bin in t2 (post_label)
dt2 = t2-Nbint2/Nbint;
Nbint3 = floor(Nbint*t3);   %%% Number of time bin in t3 (pre_label)
%Nbint1=floor(Nbint*(1-t2-t3));
Nbint_all=Nbint;
Nbint=Nbint-Nbint2-Nbint3;
% dt3 = t3-Nbint3/Nbint;
Nr = Nbint/Nbin;%%%%%%%????????%?????
dt = 1/Nbint_all;%?1??????RNA??????1/1000
dRNA = 1/Nbin;%?1??????RNA??????1/10
nRNA = 0:dRNA:Nmax;   %%% number of complete RNA molecules
Nx = length(nRNA);%mRNA????


KQ = Km*dt+diag(-kr*dt+1);%%%%%+1??%%%%?K-KINI?dt
Q0 = blkdiag0(KQ,Nx,Nstate);%%%???KQ???
% temp = cell(1,Nx);
% [temp{:}] = deal(sparse(KQ));
% Q0 = blkdiag(temp{:});

PP = zeros(Nx*Nstate,1);
PP(1:Nstate) = P0;   %%% initial state
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Propagation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t2 > 0
    ii = Nbint;
    iel = round((ii-1)/Nr);%1??bin????bin
%     if dt2 > 0
%         KQ = Km*dt2-diag(kr*dt2+1);
%         Q02 = blkdiag0(KQ,Nx,Nstate);
% %         temp = cell(1,Nx);
% %         [temp{:}] = deal(sparse(KQ));
% %         Q02 = blkdiag(temp{:});
%         QQ = Q02+sparse((iel*Nstate+1):Nx*Nstate,1:(Nx-iel)*Nstate,repmat(kr*dt2,1,Nx-iel),Nx*Nstate,Nx*Nstate);
%         PP = QQ*PP;
%     end
    
    temp = kr'*dt;%kinidt
    temp1 = temp(:,ones(1,Nx-iel));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp1 = temp1(:)';
    %Q_TEST1=sparse((iel*Nstate+1):Nx*Nstate,1:(Nx-iel)*Nstate,temp1,Nx*Nstate,Nx*Nstate);
    QQ = Q0+sparse((iel*Nstate+1):Nx*Nstate,1:(Nx-iel)*Nstate,temp1,Nx*Nstate,Nx*Nstate);
    for ii = Nbint2:-1:1
        PP = QQ*PP;
        %PP(PP<0) = 0;
    end
else
    iel = nan;
end

iel0 = iel;
for ii = Nbint:-1:1
    %ii=2
    iel = round((ii-1)/Nr);
    if iel ~= iel0
        temp = kr'*dt;
        temp1 = temp(:,ones(1,Nx-iel));
        temp1 = temp1(:)';
        QQ = Q0+sparse((iel*Nstate+1):Nx*Nstate,1:(Nx-iel)*Nstate,temp1,Nx*Nstate,Nx*Nstate);
        iel0 = iel;
    end
    PP = QQ*PP;
    %PP(PP<0) = 0;
end

%if Nbint3 > 0
%    for ii = Nbint3:-1:1
%        PP = QQ*PP;
%        PP(PP<0) = 0;
%    end
%end

%PP(PP<0) = 0;
PP = PP/sum(PP);
%tfalse = any(PP < 0);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% log likelihood calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PP0 = sum(reshape(PP,Nstate,Nx),1);
P_dtb_rna=[PP0(1),sum(reshape(PP0(2:end),[Nbin,Nmax0]),1)];
% P_dtb_rna_single=(1-alp)*[0,P_dtb_rna(2:end)];
% P_dtb_rna_nz=[0,P_dtb_rna(2:end)];P_dtb_rna_nz=P_dtb_rna_nz/sum(P_dtb_rna_nz);
% P_dtb_double=conv(P_dtb_rna_nz,P_dtb_rna_nz);
% P_dtb_rna_double=alp*P_dtb_double;
% P_dtb_rna=[P_dtb_rna(1),zeros(1,length(P_dtb_rna)-1)]+(1-P_dtb_rna(1)).*(P_dtb_rna_single+P_dtb_rna_double(1:length(P_dtb_rna_single)));
p_conv=alp*conv(P_dtb_rna,P_dtb_rna);
P_dtb_rna=(1-alp)*P_dtb_rna+p_conv(1:length(P_dtb_rna));
PP0=P_dtb_rna;
% PP0 = (PP0(round(Nbin/2):Nbin:end)+PP0(round(Nbin/2)+1:Nbin:end))/2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(pTr)
    PP0 = PP0*pTr;
end
% PP0=(1-p0)*PP0;
% PP0(1)=PP0(1)+p0;
nRNA1 = nRNA(1:Nbin:end);
% Nbin = 1;
% PP0(PP0<0) = 0;
% PP0(1) = sum(PP0(1:Nth));
% PP0(2:Nth) = 0;
PP0 = PP0/sum(PP0);

logPP = log(PP0);
% NnRNA = hist(ndata,nRNA);
%ndata=f_ob
%if tfalse
%    lld = inf;
%else
%     lld = -NnRNA*logPP';
% %     ndata(ndata > Nmax) = Nmax;
ndata = ndata(ndata <= TXmax);
%     ndata(ndata > TXmax) = TXmax;
ndata(ndata < Nth) = 0;
lld = -sum(logPP(round(ndata)+1));
    %lld=-sum(PP0(round(ndata)+1));
%end
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

