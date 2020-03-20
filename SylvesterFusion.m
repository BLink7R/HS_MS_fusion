function [ X ] = SylvesterFusion( MS,HS,b,R,subspace,Sigma,mu,Lambda_m,Lambda_h )
%SYLVESTERFUSION fuse image
%input
%   MS(data cube): MS image
%   HS(data cube): HS image
%   b: bluring kernal
%   R: spectral degrade matrix
%   subspace: the subspace X lay in
%   Sigma: a prior convariance of X
%   Mu: a prior expect of X
%   Lambda_m: the diagnol of noice covariance in MS
%   Lambda_h: the diagnol of noice covariance in HS
%output
%   X(data cube): fused image

%%
[MSrow,MScol,MSband] = size(MS);
[HSrow,HScol,HSband] = size(HS);
MSmat = img2mat(MS);
HSmat = img2mat(HS);
scale = MSrow/HSrow; % downsample factor

%% calculate a prior expect of X in subspace
% mu = imresize(HS,scale,'bicubic');
% mu = ima_interp_spline(HS,scale);
% figure(5);imshow(mu(:,:,2:4));title('mu');
invspace = pinv(subspace);
mu = invspace*img2mat(mu);

%% converse kernal b to matrix B
B = kernal2mat(b,[MSrow MScol]);
FB = fft2(B);
ratio = scale*scale;
D = reshape(FB,HSrow*HScol,ratio);  % eigen value of B
sumD = sum((D.^conj(D)),2);
Dbar = spdiags(D,HSrow*HScol*(0:-1:1-ratio),...
        sparse(MSrow*MScol,HSrow*HScol));

%% start fuse
invL_m = diag(1./Lambda_m);
invL_h = diag(1./Lambda_h);
invSig = diag(1./Sigma);

maskHS = UpsampleIMG(HS,scale);
HsBS = img2mat(ifft2(fft2(maskHS).*FB));

HLH = subspace'*invL_h*subspace;
RH = R*subspace;
C1 = HLH\(RH'*invL_m*RH+invSig);
[Q,Lambda] = eig(C1);
C3 = (HLH*Q)\(subspace'*invL_h*HsBS + RH'*invL_m*MSmat + invSig*mu);
C3 = img2mat(fft2(mat2img(C3,MSrow)));

U = zeros(size(subspace,2),MSrow*MScol);
for i = 1:size(subspace,2) % for each row
    row = C3(i,:)/Lambda(i,i);
    mid = spdiags(1./(Lambda(i,i)*ratio+sumD),0,sparse(HSrow*HScol,HSrow*HScol));
    U(i,:) = row-row*Dbar*mid*Dbar';
end

%%
X = subspace*Q*U;
X = real(ifft2(mat2img(X,MSrow)));

end

