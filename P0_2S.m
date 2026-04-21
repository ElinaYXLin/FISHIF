function P0all = P0_2S(kon,koff,kTX)

k0 = kon+koff;
k2 = kon-koff-kTX;
% k2 = -kon+koff+kTX;
k = sqrt(k2.^2+4*kon.*koff);
ka = (k+k2)/2;
kb = (k-k2)/2;
Pon = kon./k0;
Poff = koff./k0;

P0_th = (kb./k.*exp(-kb)+ka./k.*exp(ka)).*exp(-kon);
P0all = P0_th.*Pon;
P0_th = (-1./k.*exp(-kb)+1./k.*exp(ka)).*exp(-kon).*koff;
P0all = P0all+P0_th.*Pon;
P0_th = (-1./k.*exp(-kb)+1./k.*exp(ka)).*exp(-kon).*kon;
P0all = P0all+P0_th.*Poff;
P0_th = (ka./k.*exp(-ka)+kb./k.*exp(kb)).*exp(-(koff+kTX));
P0all = P0all+P0_th.*Poff;