function [ X ] = ADMMFusion( MS,HS,b,R,subspace,Xest,parameter )
%SPARSEFUSION
%input
%   MS(data cube): MS image
%   HS(data cube): HS image
%   b: bluring kernal
%   R: spectral degrade matrix
%   subspace: the subspace X lay in
%   Xest: a rough estimation of X
%   parameter: cell array of Lambda_m,Lambda_h,Lambda,Mu
%output
%   X(data cube): fused image

%%
[MSrow,MScol,MSband] = size(MS);
[HSrow,HScol,HSband] = size(HS);
MSmat = img2mat(MS);
HSmat = img2mat(HS);
scale = MSrow/HSrow; % downsample factor

[Lambda_m,Lambda_h,Lambda,Mu] = deal(parameter{:});
invL_h = diag(1./Lambda_h);
invL_m = diag(1./Lambda_m);

%% calculate a prior mean of X in subspace
subspace_dim = size(subspace,2);
invspace = pinv(subspace);
Uest = invspace*img2mat(Xest);
Uest = mat2img(Uest,MSrow);

%% converse kernal b to matrix B
B = kernal2mat(b,[MSrow MScol]);
FB = fft2(B); % eigen values of B
IB = 1./(2+FB.*conj(FB)); % eigen values of (B*B^T+2I)^-1

%% iteration
fprintf('start ADMM step\n');
V1 = BlurIMG(Uest,b);
V2 = img2mat(Uest);
V3 = V2;
G1 = zeros(size(V1));
G2 = zeros(size(V2));
G3 = G2;
for ADMM_ite=1:30
    fprintf('ADMM %d iteration\n',ADMM_ite);
    fprintf('Start SALSA sub-iteration\n');
    for SALSA_ite=1:5
        fprintf('SALSA %d iteration\n',SALSA_ite);
        Uopt = ifft2(fft2(V1+G1).*conj(FB))+mat2img(V2+G2+V3+G3,MSrow);
        Uopt = real(ifft2(fft2(Uopt).*IB));
        Umat = img2mat(Uopt);
        % V1
        v1 = BlurIMG(Uopt,b)-G1;
        V1 = v1;
        V1mat = (subspace'*invL_h*subspace+Mu*eye(HSband))...
            \(subspace'*invL_h*HSmat+Mu*img2mat(DownsampleIMG(v1,scale)));
        V1(1:scale:MSrow,1:scale:MScol,:)=mat2img(V1mat,HSrow);
        %     figure(6);imshow(V1(:,:,2:4));title('V1');
        % V2
        v2 = Umat-G2;
        HR = R*subspace;
        V2 = (HR'*invL_m*HR+Mu*eye(subspace_dim))\(HR'*invL_m*MSmat+Mu*v2);
        % V3
        v3 = Umat-G3;
        V3 = (Lambda+Mu)\(Lambda*img2mat(Uest)+Mu*v3);
        % G
        G1 = -v1+V1;
        G2 = -v2+V2;
        G3 = -v3+V3;
    end
    % print result
    X = real(subspace*Umat);
    X = mat2img(X,MSrow);
    figure(5);imshow(X(:,:,2:4));title('X');
    Uest = Uopt;
end
fprintf('ADMM ok\n');
%%
X = real(subspace*Umat);
X = mat2img(X,MSrow);