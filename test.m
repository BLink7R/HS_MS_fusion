% clc,clear,close all
% load('matlab.mat')
% origin = double(imread('duck.jpg'));
% origin = origin(1:600,1:600,:);
% origin = origin/max(origin(:));
% 
% kernal = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1];
% kernal = kernal/sum(kernal(:));
% HS = BlurIMG(origin,kernal);
% HS = DownsampleIMG(HS,2);
% % HSim = imnoise(HSim,'gaussian',0.01);
% MS = mat2img([0.4 0.3 0.3]*img2mat(origin),size(origin,1));
% % MSim = imnoise(MSim,'gaussian',0.01);
% figure(1);imshow(origin);title('origin');

%%
addpath('./gadget');
addpath('./dicLearn');
load('data.mat');
% HS = MS;
% MS = PAN;
HS = MS(251:400,151:300,:);
MS = PAN(1001:1600,601:1200,:);

[HSrow,HScol,HSband]=size(HS);
[MSrow,MScol,MSband]=size(MS);
scale=MSrow/HSrow;

figure(1);imshow(HS(:,:,2:4));title('HS');
figure(2);imshow(MS);title('MS');

%% parameter
Lambda_b = 1e5; % regularization parameter
Lambda_r = 1e1; % regularization parameter
rB       = 7; % size of B

subspace_dim = 4;

Lambda_m = 10*ones(1,MSband);
Lambda_h = 10*ones(1,HSband);

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
[V,~] = svd(img2mat(HS));
subspace = V(:,1:subspace_dim);
% subspace = eye(4);
% [subspace, ~, ~, ~]=idHSsub(HS,'PCA',1,dimension);

%% rough estimation of X
% Xest = imresize(HS,scale,'bicubic');
Xest = ima_interp_spline(HS,scale);

%% image fusion-SylvesterFusion
% Sigma = 1e3*ones(1,subspace_dim);
% X = SylvesterFusion(MS,HS,B,R,subspace,Sigma,Xest,Lambda_m,Lambda_h);
% figure(5);imshow(X(:,:,2:4));title('X');

%% image fusion-SparseFusion
Lambda = 0.1;
Mu = 0.001;
patch_size = [6 6];
basis_num = [16 16];
parameter={Lambda_m,Lambda_h,Lambda,Mu};
X = SparseFusion(MS,HS,B,R,subspace,Xest,parameter,patch_size,basis_num);
figure(5);imshow(X(:,:,2:4));title('X');

%% check result
[snr1,snr2]=CheckResult(HS,MS,X,R,B);
fprintf('snr in HS is %12.5f,\nsnr in MS is %12.7f\n',snr1,snr2);