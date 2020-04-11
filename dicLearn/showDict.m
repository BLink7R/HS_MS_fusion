function showDict( D,lin,col )
%SHOWDICT 
%input
%   D: the dict
[basis_size,basis_num] = size(D);
basis_size = round(sqrt(basis_size));
D = D-min(D(:));
D = D./max(D(:));
D = reshape(D,basis_size,basis_size,basis_num);
% draw on a board, with a border of 2px
board = ones((basis_size+2)*col-2,(basis_size+2)*lin-2);
for i=0:lin-1
    for j=0:col-1
        up    = i*(basis_size+2)+1;
        down  = up+basis_size-1;
        left  = j*(basis_size+2)+1;
        right = left+basis_size-1;
        board(up:down,left:right) = D(:,:,i*col+j+1);
    end
end
imshow(board);