function [ B,R ] = EstimateBR( MS,HS,band_map,HS_bandID,Lambda_b,Lambda_r,rB )
%ESTIMATEBR estimate degrade matrix B and R
%input
%   MS(data cube): MS image
%   HS(date cube): HS image
%       (the row number and colum number of HS should be divisors of MS's)
%   band_map(cell): the band numbers of HS image corresponding to each MS band
%   HS_bandID: record whether HS bands are contiguous
%   Lambda_b,Lambda_h: regularization parameter
%   sizeB: decide the size of B (2*sizeB+1,2*sizeB+1)
%output
%   B,R: degrade matrix

%%
[MSrow,MScol,MSband] = size(MS);
[HSrow,HScol,HSband] = size(HS);
scale = MSrow/HSrow; % downsample factor

%% blur HS and MS using a strong bluring kernal
r = 1; % parameter
bluredMS = BlurIMG(MS,fspecial('average',[2*r*scale+1,2*r*scale+1]));
bluredHS = BlurIMG(HS,fspecial('average',[2*r+1,2*r+1]));
bluredMS = img2mat(DownsampleIMG(bluredMS,scale));
bluredHS = img2mat(bluredHS);

%% estimate R
R = zeros(MSband,HSband);
for i=1:MSband
    bands_list = band_map{i};
    bands_l    = length(bands_list);
    bands      = bluredHS(bands_list,:);
    contiguous = (diff(HS_bandID(bands_list))==1);
    D = sparse(bands_l,bands_l-1);
    if bands_l>1, D = spdiags([-contiguous contiguous],[0 -1],D); end
    R(i,bands_list) = bluredMS(i,:)*bands'/(bands*bands'+Lambda_r*(D*D'));
end

%% estimate B
sizeB = (2*rB+1)^2;
padMS = padarray(MS,[rB rB],'circular'); % circulant boundary condition
% padMS = padarray(MS,[sizeB sizeB]); %zero boundary condition
% degraded HS
downHS = mat2img(R*img2mat(HS),HSrow);
% insert 0 in HS to match MS
maskHS = UpsampleIMG(downHS,scale);
% diff matrixes
Dh = spdiags([ones(sizeB,1) -ones(sizeB,1)],[0 2*rB+1],sparse(sizeB-(2*rB+1),sizeB));
Dv = kron(eye(2*rB+1),...
    spdiags([ones(2*rB,1) -ones(2*rB,1)],[0 1],sparse(2*rB,2*rB+1)));

sum1 = Lambda_b*(Dh'*Dh+Dv'*Dv);
sum2 = zeros(sizeB,1);
% calculate each pixel as a vector of 1*MSband
for y=rB+1:scale:MSrow+rB
    for x=rB+1:scale:MScol+rB
        block = img2mat(padMS(y-rB:y+rB,x-rB:x+rB,:)); % pixels affected by bluring kernal
        sum1  = sum1+block'*block;
        sum2  = sum2+block'*img2mat(maskHS(y-rB,x-rB,:));
    end
end
b = sum1\sum2;

%% normalize
sumB = sum(b(:));
b = b/sumB;
R = R/sumB;
B = rot90(reshape(b,2*rB+1,2*rB+1),2);

end