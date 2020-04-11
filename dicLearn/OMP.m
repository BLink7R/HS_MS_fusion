function A = OMP( D,patch,L )
%OMP orthogonal matching pursuit
%input
%   D: the dict
%   patch: patch_size*patch_num matrix
%   L: sparse parameter
%output
%   A: sparse code

%% init
basis_num = size(D,2);
patch_num = size(patch,2);
A = sparse(basis_num,patch_num);

%% work
% use=zeros(1,basis_num);
% SNR=zeros(1,patch_num);
for i=1:patch_num % for each patch
    p = patch(:,i);
    R = p; % residual
    support = zeros(1,L);
%     reshuffle = randperm(basis_num);
    for j=1:L
        proj = abs(D'*R);
        [~,support(j)] = max(proj);
        % update support
%         Alpha = 0.6;
%         for k=1:basis_num
%             ind=reshuffle(k);
%             if proj(ind)>Alpha*maxi
%                 support(j)=ind;
%                 break;
%             end
%         end
        Dsup = D(:,support(1:j));
        % update sparse code
        a=pinv(Dsup)*p;
        R = p-Dsup*a;
    end
%     SNR(i)=snr(p,R);
    A(support(1:j),i) = a;
end
% disp(max(use));

