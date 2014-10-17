fname = 'data\IMG_4305';

r = imread([fname,'_R_corr.pgm']);
g = imread([fname,'_G_corr.pgm']);
b = imread([fname,'_B_corr.pgm']);

rgb = zeros(size(r,1), size(r,2), 3);
rgb = uint8(rgb);
rgb(:,:,1) = r;
rgb(:,:,2) = g;
rgb(:,:,3) = b;

imwrite(rgb, [fname, '_RGB_corr.bmp']);

%i = imread([fname, '_RGB.bmp']);

%figure; imshow(i);
figure; imshow(rgb);