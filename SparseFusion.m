function [ X ] = SparseFusion( MS,HS,b,R,subspace,Xest,parameter,patch_size,basis_num )
%SPARSEFUSION
%input
%   MS(data cube): MS image
%   HS(data cube): HS image
%   b: bluring kernal
%   R: spectral degrade matrix
%   subspace: the subspace X lay in
%   Xest: a rough estimation of X
%   parameter: cell array of Lambda_m,Lambda_h,Lambda,Mu
%   patch_size
%   basis_num: sqrt of basis number
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

%% calculate a prior expect of X in subspace
subspace_dim = size(subspace,2);
invspace = pinv(subspace);
Uest = invspace*img2mat(Xest);
Uest = mat2img(Uest,MSrow);

%% dict learning
patch_num = MSrow*MScol/prod(patch_size);
% D = zeros(prod(patch_size),prod(basis_num),HSband);
% A = cell(HSband,1);
load('./DA6_6.mat');
patch = zeros(prod(patch_size),patch_num,HSband);
for i=1:HSband
%     fprintf('start dictionary learning using KSVD in %d band\n',i);
%     ind = 1;
%     for y=1:patch_size(1):MSrow
%         for x=1:patch_size(2):MScol
%             p = Uest(y:y+patch_size(1)-1,x:x+patch_size(2)-1,i);
%             patch(:,ind,i) = reshape(p,prod(patch_size),1);
%             ind = ind+1;
%         end
%     end
%     [D(:,:,i),A{i}] = K_SVD(patch(:,:,i),6,basis_num(1));
%     fprintf('KSVD OK\n');
    patch(:,:,i) = D(:,:,i)*A{i};
    ind = 1;
    for y=1:patch_size(1):MSrow
        for x=1:patch_size(2):MScol
            Uest(y:y+patch_size(1)-1,x:x+patch_size(2)-1,i)=reshape(patch(:,ind,i),patch_size);
            ind = ind+1;
        end
    end
end
Xest = subspace*img2mat(Uest);
Xest = mat2img(Xest,MSrow);
figure(4); imshow(Xest(:,:,[2,3,4])); title('rough estimation');

%% converse kernal b to matrix B
B = kernal2mat(b,[MSrow MScol]);
FB = fft2(B); % eigen values of B
IB = 1./(2+FB.*conj(FB)); % eigen values of (B*B^T+2I)^-1

%% iteration
V1 = BlurIMG(Uest,b);
V2 = img2mat(Uest);
V3 = V2;
G1 = zeros(size(V1));
G2 = zeros(size(V2));
G3 = G2;
for i=1:10
    Uopt = ifft2(fft2(V1+G1).*conj(FB))+mat2img(V2+G2+V3+G3,MSrow);
    Uopt = real(ifft2(fft2(Uopt).*IB));
    % use dict learning to generate a better Uopt
    for band=1:HSband
        ind = 1;
        for y=1:patch_size(1):MSrow
            for x=1:patch_size(2):MScol
                support = A{band}(:,ind)~=0;
                d = D(:,support);
                block = Uopt(y:y+patch_size(1)-1,x:x+patch_size(2)-1,band);
                block = d*pinv(d)*reshape(block,prod(patch_size),1);
                Uopt(y:y+patch_size(1)-1,x:x+patch_size(2)-1,band) = reshape(block,patch_size);
                ind = ind+1;
            end
        end
    end
    Umat = img2mat(Uopt);
    % print result
    X = real(subspace*Umat);
    X = mat2img(X,MSrow);
    figure(5);imshow(X(:,:,2:4));title('X');
    % V1
    v1 = BlurIMG(Uopt,b)-G1;
    V1 = v1;
    V1mat = (subspace'*invL_h*subspace+Mu*eye(HSband))...
        \(subspace'*invL_h*HSmat+img2mat(DownsampleIMG(v1,scale)));
    V1(1:scale:MSrow,1:scale:MSrow,:)=mat2img(V1mat,HSrow);
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
    
    Uest = Uopt;
end

%%
X = real(subspace*img2mat(Uopt));
X = mat2img(X,MSrow);

