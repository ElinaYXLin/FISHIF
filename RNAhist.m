function y = RNAhist(beta,x)

A1 = beta(1);
A2 = beta(2);
x02 = beta(3);
sigma2 = beta(4);
A3 = beta(5);
x03 = beta(6);
sigma3 = beta(7);

y = A1*x.^(-2)+A2*exp(-(x-x02).^2./2./sigma2.^2)+A3*exp(-(x-x03).^2./2./sigma3.^2);
