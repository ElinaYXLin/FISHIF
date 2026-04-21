% fmat=moviein(1);
% fmat(:,j)=getframe;
function GifMaker(fmat,save_path)
% save_path='Y:\meeting\';
for j=1:size(fmat,2)
    im=frame2im(fmat(:,j));
    [I,map]=rgb2ind(im,256);
    if j==1
        imwrite(I,map,[save_path,'out.gif'],'gif','loopcount',inf,'DelayTime',0.2)
    else
        imwrite(I,map,[save_path,'out.gif'],'gif','writemode','append','DelayTime',0.2)
    end
end