%% load data
addpath('./gadget');
load('data.mat')
HS = MS;
MS = PAN;

[HSrow,HScol,HSband]=size(HS);
[MSrow,MScol,MSband]=size(MS);
scale = MSrow/HSrow;

matHS = img2mat(HS);
matMS = img2mat(MS);

maskHS = UpsampleIMG(HS,scale);
figure(1);imshow(1.5*HS(:,:,2:4));title('HS');
figure(2);imshow(1.5*MS);title('MS');

%% parameter
rB = 8;
sizeB = (2*rB+1)^2;
subspace_dim = 4;

Lambda_b = 1e6; % regularization parameter of B
Lambda_b2 = 1e1;
Lambda_r = 1e3; % regularization parameter of R


Lambda_m = 1*ones(1,MSband); % covariance of noice in MS
Lambda_h = 1*ones(1,HSband); % covariance of noice in HS
Sigma = 1e3*ones(1,subspace_dim);

%% band map
band_map=cell(1,MSband);
HS_bandID = [1 2 3 4]';
[~,band_map{1}] = intersect(HS_bandID, 1:4);

%% diff matrix
% for each row of R
Dr = cell(1,MSband);
for i=1:MSband
    bands_list = band_map{i};
    bands_l    = length(bands_list);
    contiguous = (diff(HS_bandID(bands_list))==1);
    Dr{i} = sparse(bands_l,bands_l-1);
    if bands_l>1, Dr{i} = spdiags([1 -1].*contiguous,[0 -1],Dr{i}); end
    Dr{i} = Lambda_r*Dr{i}*Dr{i}';
end

% for b
Dh = spdiags([1 -1].*ones(sizeB,1),[0 2*rB+1],sparse(sizeB-(2*rB+1),sizeB));
Dv = kron(eye(2*rB+1),...
    spdiags([1 -1].*ones(2*rB,1),[0 1],sparse(2*rB,2*rB+1)));
Dl = 2*eye(sizeB)-ones(sizeB);
Db = Lambda_b*(Dh'*Dh+Dv'*Dv+Lambda_b2*Dl);
%% subspace identification
% [V,~] = svd(matHS);
% subspace = V(:,1:subspace_dim);
[subspace, ~, ~, ~]=idHSsub(HS,'PCA',1,subspace_dim);
% subspace = eye(subspace_dim);

%% iteration initial value
Xk = ima_interp_spline(HS,scale);
X0 = Xk;
figure(3);imshow(1.5*Xk(:,:,[2 3 4]));title('X0');
% bk = fspecial('average',[2*rB+1 2*rB+1]);
% Rk = fspecial('average',[MSband HSband]);

%% iteration
ITE_TIME = 0;
while ITE_TIME<10
    ITE_TIME = ITE_TIME+1;
    matXk = img2mat(Xk);
    % solve R
    Rk = zeros(MSband,HSband);
    for i=1:MSband
        bands_list = band_map{i};
        bands = matXk(bands_list,:);
        Rk(i,bands_list) = matMS(i,:)*bands'/(bands*bands'+Dr{i});
    end
    disp('r ok');
    % solve b (kernal of B)
    sum1 = Db;
    sum2 = -Lambda_b2*Lambda_b*ones(sizeB,1);
    padX = padarray(Xk,[rB rB],'circular');
    for y=1:scale:MSrow
        for x=1:scale:MScol
            block = img2mat(padX(y:y+2*rB,x:x+2*rB,:)); % pixels affected by bluring kernal
            sum1  = sum1+block'*block;
            sum2  = sum2+block'*img2mat(maskHS(y,x,:));
        end
    end
    % show bluring kernal
    bk = rot90(reshape(sum1\sum2,2*rB+1,2*rB+1));
    disp('b ok');
    figure(5);imagesc(bk);axis image;axis off;
    set(gca,'FontSize',15);
    colorbar;title('Estimated Spatial Blurring');
    % solve X
    Xk = SylvesterFusion(MS,HS,bk,Rk,subspace,Sigma,Xk,Lambda_m,Lambda_h);
%     Sigma = Sigma.*0.5;
    disp('X ok');
    figure(4);imshow(1.5*Xk(:,:,2:4));title('X');
    [snr1,snr2]=CheckResult(HS,MS,Xk,Rk,bk);
    fprintf('snr in HS is %12.5f,\nsnr in MS is %12.7f\n',snr1,snr2);
end
