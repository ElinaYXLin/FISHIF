function errorbarxy(x,y,lx,ly,ux,uy,linecol,errorcol,varargin)
%This function allows the user to plot the graph of x against y, along with both x and y errorbars.
%For the x and y errors it is possible to input both lower (lx and ly)  and upper  (ux and uy) values for the
%errors at a particular point.  If the upper values are not specified then the program assumes the errors 
%are symmetrical and use the lower values.  it is also possible to specify the plot line colour, marker, and 
%linestyle using the standard 'plot' command notation in the input variable 'linecol'.  Also the line colour for 
%the errobars can be specified in the variable 'errorcol'.  It is important to note that if these colour options 
%are to be used and any of the error limit vectors are empty then they should not be excluded, but presented 
%in a [] form signifying an empty vector.
%
%James Rooney,  17 October 2003
temp_hold = ishold;

if exist('linecol','var')==0 | isempty(linecol)
    linecol='b';
end


if exist('errorcol','var')==0 | isempty(errorcol)
    errorcol='r';
end

%axe = axes;
hold on

xw=(max(x)-min(x))/100;
yw=(max(y)-min(y))/100;


lye=exist('ly','var');
lxe=exist('lx','var');
uye=exist('uy','var');
uxe=exist('ux','var');

if lye+lxe+uye+uxe==0 | isempty(lx) & isempty(ux) & isempty(ly) & isempty(uy)
    return
end

if uye==0 | isempty(uy)
    uy=ly;
end

if uxe==0 | isempty(ux)
    ux=lx;
end

ex_x = zeros(0);
ex_y = zeros(0);
ey_x = zeros(0);
ey_y = zeros(0);

for t=1:length(x)
    if ~isempty(ux)
    %x errorbars
        ex_x = cat(1,ex_x,[x(t)-lx(t);x(t)+ux(t)],nan,[x(t)-lx(t);x(t)-lx(t)],nan,[x(t)+ux(t);x(t)+ux(t)],nan);
        ex_y = cat(1,ex_y,[y(t);y(t)],nan,[y(t)-yw;y(t)+yw],nan,[y(t)-yw;y(t)+yw],nan);
 
    end
    
    if ~isempty(uy)
    %y errorbars
        ey_x = cat(1,ey_x,[x(t);x(t)],nan,[x(t)-xw;x(t)+xw],nan,[x(t)-xw;x(t)+xw],nan);
        ey_y = cat(1,ey_y,[y(t)-ly(t);y(t)+uy(t)],nan,[y(t)-ly(t);y(t)-ly(t)],nan,[y(t)+uy(t);y(t)+uy(t)],nan);

    end    
end
e_x = cat(1,ex_x,ey_x);
e_y = cat(1,ex_y,ey_y);
%line(e_x,e_y,'color',errorcol);
if ~isempty(varargin) && any(strcmp(varargin{1},{'+','o','*','.','x','s','d','^','v','>','<','p','h'}))
    special_mark = varargin{1};
    special_line = 'none';
elseif ~isempty(varargin) && any(strcmp(varargin{1},{'-','--',':','-.'}))
    special_mark = 'none';
    special_line = varargin{1};
else
    special_mark = 'none';
    special_line = '-';
end
h = errorbar(x,y,zeros(size(x)),'color',linecol,varargin{2:end});
h2 = get(h,'children');
if ~isempty(h2)
    set(h2(1),'Marker','none','LineStyle','none');
    set(h2(2),'Marker','none','LineStyle','none');
end
if ~isempty(varargin) && length(varargin) > 1
    try
        plot(e_x,e_y,'color',errorcol,'parent',h);
    %     plot(e_x,e_y,'color',errorcol,'parent',h,varargin{2:end});
    catch
        plot(e_x,e_y,'color',errorcol);
    end
else
    try
        plot(e_x,e_y,'color',errorcol,'parent',h);
    catch
        plot(e_x,e_y,'color',errorcol);
    end
end

try
    plot(x,y,'Marker',special_mark,'color',linecol,'LineStyle',special_line,'parent',h);
catch
    plot(x,y,'Marker',special_mark,'color',linecol,'LineStyle',special_line);
end



% hg = hggroup;
% set(h,'Parent',hg)
% set(get(get(hg,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% set(hg,'HandleVisibility','on')
% set(hg,'DisplayName','111')

% hg = hgtransform;%('Parent',ax);
% set(h,'Parent',hg)

if ~temp_hold
    hold off
end
    
