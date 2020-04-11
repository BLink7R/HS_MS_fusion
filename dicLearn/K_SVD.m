function [D,A] = K_SVD(patch,K,basis_num)
%K-SVD
%input
%   patch: patch_size*patch_num matrix
%   K: sparse parameter
%   basis_num: sqrt of basis number
%output
%   D: optimized dictionary
%   A: optimized sparse code

%%
basis_size = round(sqrt(size(patch,1)));
D=overcompleteDCT(basis_size,basis_num);
% D=rand(basis_size^2,basis_num^2);
% normlize to 1
for i=2:basis_num^2, D(:,i)=D(:,i)/norm(D(:,i),2); end

%% work
for iter=1:15
    fprintf('KSVD %d iteration\n',iter);
%     figure(2);showDict(D,basis_num,basis_num);
    % orthogonal matching pursuit
    A=OMP(D,patch,K);
    for i=1:basis_num^2
        % get support
        support=find(A(i,:));
        if isempty(support), continue; end
        % calculate residual
        A(i,:)=0;
        Ei=patch(:,support)-D*A(:,support);
        % use svd to find rank-1 approx
        [U,L,V]=svds(Ei,1);
        A(i,support)=L*V*norm(U,2);
        D(:,i)=U/norm(U,2);
    end
end
