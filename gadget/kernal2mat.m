function mat = kernal2mat( kernal,mat_size )
%KERNAL2MAT converse a bluring kernal to matrix
mat = zeros(mat_size(1),mat_size(2));
ry = floor((size(kernal,1)+1)/2);
rx = floor((size(kernal,2)+1)/2);
midy = round((mat_size(1)+1)/2);
midx = round((mat_size(2)+1)/2);
mat(midy-ry+1:midy-ry+size(kernal,1),midx-rx+1:midx-rx+size(kernal,2)) = kernal;
mat = ifftshift(mat);
end

