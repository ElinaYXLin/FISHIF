function n2D = hist3w(n_raw,bin2D,w2D)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to calculation 2D histogram with adjustable bin width %%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% n2D (output): 2D histogram result;
%% n_raw (input): raw data (value in i, value in j);
%% bin2D (input): 2D bin {bin in i, bin in j};
%% w2D (input): 2D bin width (w in i, w in j);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin1 = bin2D{1};
bin2 = bin2D{2};
N1 = length(bin1);
N2 = length(bin2);
w1 = w2D(1)/2;
w2 = w2D(2)/2;
n2D = zeros(N1,N2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2D Histogram generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:N1
    for jj = 1:N2
        n2D(ii,jj) = nnz((n_raw(:,1) >= (bin1(ii)-w1)) & (n_raw(:,1) <= (bin1(ii)+w1)) & (n_raw(:,2) >= (bin2(jj)-w2)) & (n_raw(:,2) <= (bin2(jj)+w2)));
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


