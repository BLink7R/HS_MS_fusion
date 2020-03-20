function img = mat2img( mat,row )
%MAT2IMG 
img = reshape(mat',row,size(mat,2)/row,size(mat,1));
end

