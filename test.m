%%
addpath('./gadget');
addpath('./dicLearn');
load('data.mat');
HS = MS;
MS = PAN;
% HS = MS(251:400,151:300,:);
% MS = PAN(1001:1600,601:1200,:);

[HSrow,HScol,HSband]=size(HS);
[MSrow,MScol,MSband]=size(MS);
scale=MSrow/HSrow;

figure(1);imshow(HS(:,:,2:4));title('HS');
figure(2);imshow(MS);title('MS');

%% parameter
Lambda_b = 1; % regularization parameter
Lambda_r = 1; % regularization parameter
rB       = 7; % size of B

subspace_dim = 4;

%% estimate B and R
band_map=cell(1,MSband);
non_del_bands = [1 2 3 4];
[~,band_map{1}] = intersect(non_del_bands, 1:4);
fprintf('start estimating degrade matrix B and R...\n');
[B,R] = EstimateBR(MS,HS,band_map,non_del_bands',Lambda_b,Lambda_r,rB);
figure(3);imagesc(B);axis image;axis off;
set(gca,'FontSize',15);
colorbar;title('Estimated Spatial Blurring');
fprintf('estimate OK\n');

%% subspace identification
% [V,~] = svd(img2mat(HS));
% subspace = V(:,1:subspace_dim);
subspace = eye(4);
% [subspace, ~, ~, ~]=idHSsub(HS,'PCA',1,dimension);

%% rough estimation of X
% Xest = imresize(HS,scale,'bicubic');
Xest = ima_interp_spline(HS,scale);

%% image fusion-SylvesterFusion
Lambda_m = 1*ones(1,MSband);
Lambda_h = 1000*ones(1,HSband);
Sigma = 1e2*ones(1,subspace_dim);
X = SylvesterFusion(MS,HS,B,R,subspace,Sigma,Xest,Lambda_m,Lambda_h);

%% image fusion-ADMMFusion
Lambda_m = 1*ones(1,MSband);
Lambda_h = 1*ones(1,HSband);
Lambda = 0.1;
Mu = 0.001; % ADMM step length
parameter = {Lambda_m,Lambda_h,Lambda,Mu};
X = ADMMFusion(MS,HS,B,R,subspace,Xest,parameter);

%% check result
figure(5);imshow(X(:,:,2:4));title('X');
[snr1,snr2]=CheckResult(HS,MS,X,R,B);
fprintf('snr in HS is %12.5f,\nsnr in MS is %12.7f\n',snr1,snr2);