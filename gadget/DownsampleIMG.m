function IMG = DownsampleIMG( IMG,ratio )
%UPSAMPLEIMG 
IMG = permute(downsample(permute(downsample(IMG,ratio),[2 1 3]),ratio),[2 1 3]);
end

