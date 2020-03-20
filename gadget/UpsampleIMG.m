function IMG = UpsampleIMG( IMG,ratio )
%UPSAMPLEIMG 
IMG = permute(upsample(permute(upsample(IMG,ratio),[2 1 3]),ratio),[2 1 3]);
end

