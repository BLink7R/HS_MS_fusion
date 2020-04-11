function DCT = overcompleteDCT(basis_size,basis_num)
%overcompleteDCT
%input
%   basis_size: each basis is basis_size*basis_size
%   basis_num: generate basis_num*basis_num basis
%output
%   DCT: overcomplete DCT dict
V = sqrt(2 / basis_size)*cos((0:basis_size-1)'*(0:basis_num-1)*pi/basis_size/2); 
V(1,:) = V(1,:) / sqrt(2);
DCT=kron(V,V);