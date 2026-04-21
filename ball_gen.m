function se = ball_gen(vtpsize)

se0 = strel('ball',vtpsize,vtpsize);
se = false((2*vtpsize+1),(2*vtpsize+1),(2*vtpsize+1));
se(:,:,vtpsize+1) = getnhood(se0);
he3 = getnhood(se0);
for kk = 1:(2*vtpsize+1)
    se(:,:,kk) = se(:,:,vtpsize+1).*(abs(he3-kk) > 0);
end