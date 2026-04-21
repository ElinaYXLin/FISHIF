clear all

h0 = get(gcf,'Children');
for ii = 1:length(h0)
    h1 = get(h0(length(h0)-ii+1),'Children');
    eval(['xx',num2str(ii,'%02u'),' = get(h1(3),''XData''); ']);
    eval(['yy',num2str(ii,'%02u'),' = get(h1(3),''YData''); ']);
%     axes(h0(length(h0)-ii+1))
%     ylim([0,40])
end
    