function Gxy = Con_calculator(Nx,Ny,sigma_x,sigma_y,varargin)

% % Nx = 10;
% % Ny = 10;
% % sigma_x = 1;
% % sigma_y = 1.37;
if isempty(varargin)
    Gx = 0;
    Gy = 0;
    for Id = 0:(Nx-1)
        Gx = Gx + (Nx-Id)*exp(-Id.^2./4./(sigma_x.^2));
    end
    for Id = 0:(Ny-1)
        Gy = Gy + (Ny-Id)*exp(-Id.^2./4./(sigma_y.^2));
    end
    Gxy = Gx*Gy./Nx./Ny;
else
    work_mask = double(varargin{1} ~= 0);
    work_conv = conv2(work_mask,work_mask);
    work_num = work_conv(ceil(size(work_conv,1)./2):end,ceil(size(work_conv,2)./2):end);

    [YY,XX] = meshgrid(0:(size(work_mask,2)-1),0:(size(work_mask,1)-1));

    Gxy = sum(sum(work_num.*exp(-XX.^2./4./(sigma_x.^2)-YY.^2./4./(sigma_y.^2))))./sum(sum(work_mask));
end
