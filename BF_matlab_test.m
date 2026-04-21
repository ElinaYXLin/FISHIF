% BF-MATLAB
data = bfopen([input_folder,input_name]);
seriesCount = size(data, 1);
series1 = data{1, 1};
series2 = data{2, 1};
series3 = data{3, 1};
metadataList = data{1, 2};
% etc
series1_planeCount = size(series1, 1);
series1_plane1 = series1{1, 1};
series1_label1 = series1{1, 2};
series1_plane2 = series1{2, 1};
series1_label2 = series1{2, 2};
series1_plane3 = series1{3, 1};
series1_label3 = series1{3, 2};

series1_colorMap1 = data{1, 3}{1, 1};
figure('Name', series1_label1);
if isempty(series1_colorMap1)
  colormap(gray);
else
  colormap(series1_colorMap1);
end
imagesc(series1_plane1);


cmap = gray(256);
for p = 1 : size(series1, 1)
  M(p) = im2frame(uint8(series1{p, 1}), cmap);
end
if feature('ShowFigureWindows')
  movie(M);
end