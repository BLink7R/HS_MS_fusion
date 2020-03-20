function mat = img2mat( img )
%IMG2MAT 
mat = reshape(img,size(img,1)*size(img,2),size(img,3))';
end

