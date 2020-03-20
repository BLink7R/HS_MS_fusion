function [ blured ] = BlurIMG( IMG,kernal )
%BLURIMG blur a image
% blured = conv2(IMG,kernal,'same');
blured = real(ifft2(fft2(kernal2mat(kernal,size(IMG))).*fft2(IMG)));
end
