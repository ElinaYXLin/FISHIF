function tiffwrite0(A,filename)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to write multi-channel TIFF file %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(A,3) <= 4
    imwrite(A,filename,'Compression','none')
else
    t = Tiff(filename,'w');   % Create a new image.
    t.setTag('Photometric',Tiff.Photometric.Separated);   % Set the photometric interpretation to be CMYK.

    t.setTag('ImageWidth',size(A,2));   % Set dimension X of the image in the file.
    t.setTag('ImageLength',size(A,1));   % Set dimension Y of the image in the file.
    
    info0 = whos('A');
    t.setTag('BitsPerSample',info0.bytes/numel(A)*8);   % Set the Bits of the image in the file.
    
    t.setTag('SamplesPerPixel',size(A,3));   % Set dimension Z of the image in the file.

    t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);   %%% Set the ExtraSample mode
    
    t.write(A);
    t.close();

end