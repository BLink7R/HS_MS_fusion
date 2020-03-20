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
load('data.mat')
HS = MS;
MS = PAN;

[HSrow,HScol,HSband]=size(HS);
[MSrow,MScol,MSband]=size(MS);
scale=MSrow/HSrow;

% HS = NormColor(HS);
% MS = NormColor(MS);
% HS = HS/max(HS(:));
% MS = MS/max(MS(:));

figure(2);imshow(HS(:,:,2:4));title('HS');
figure(3);imshow(MS);title('MS');

%% estimate B and R
Lambda_b = 1e5; % regularization parameter
Lambda_r = 1e1; % regularization parameter
rB       = 7; % size of B

band_map=cell(1,MSband);
non_del_bands = [1 2 3 4];
[~,band_map{1}] = intersect(non_del_bands, 1:4);

[B,R] = EstimateBR(MS,HS,band_map,non_del_bands',Lambda_b,Lambda_r,rB);

figure(5);imagesc(B);axis image;axis off;
set(gca,'FontSize',15);
colorbar;title('Estimated Spatial Blurring');
%% subspace identification
dimension = 4; % parameter
% [V,~] = svd(img2mat(HS));
% subspace = V(:,1:dimension);
[subspace, ~, ~, ~]=idHSsub(HS,'PCA',1,dimension);
% subspace = eye(4);

%% image fusion
Lambda_m = 10*ones(1,MSband);
Lambda_h = 10*ones(1,HSband);
Sigma = 1e3*ones(1,dimension);
mu = ima_interp_spline(HS,scale);
X = SylvesterFusion(MS,HS,B,R,subspace,Sigma,mu,Lambda_m,Lambda_h);
figure(4);imshow(X(:,:,2:4));title('X');

%% check result
[snr1,snr2]=CheckResult(HS,MS,X,R,B);
fprintf('snr in HS is %12.5f,\nsnr in MS is %12.7f\n',snr1,snr2);