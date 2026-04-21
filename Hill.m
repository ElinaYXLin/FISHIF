function y = Hill(beta,x)
 y = beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);
