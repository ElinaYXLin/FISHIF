function bw_out = bwthicken(bw_in,nr)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to thicken the binary image without connecting objects %%%%%%
%% It is an improvement for bwmorph('thicken'), which has bugs. %%%%%%%%%%%
%% bw_out: Output binary image.
%% bw_in: Input binary image.
%% nr: thicken times.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bw_in = logical(bw_in);
bw_out = bwmorph(bw_in,'thicken',nr);

if max(max(bwlabel(bw_in))) ~= max(max(bwlabel(bw_out)))
    D= bwdist(~bw_out);
    g2 = imimposemin(-D,bw_in);
    bw_out = ~(watershed(g2) == 0) & bw_out;
end


