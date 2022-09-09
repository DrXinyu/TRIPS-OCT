function M = H2M(H)
% converts and input matrix H into its corresponding Mueller-matrix,
% following Cloude.

N = size(H,3);

A = 1/sqrt(2)*[1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0];

for tot = 1:N
    H([2,4,6,8,9,11,13,15]) = H([9,11,13,15,2,4,6,8]);
    M = A*H(:,:,tot)*A';
end
    

% N = size(H,3);
% 
% sig(:,:,1) = eye(2);
% sig(:,:,2) = diag([1,-1]);
% sig(:,:,3) = [0 1;1 0];
% sig(:,:,4) = [0 -1i; 1i 0];
% 
% M = zeros(size(H));
% 
% for tot = 1:N
%     for ind = 1:4
%         for jnd = 1:4
%             M(ind,jnd,tot) = 1/2*trace(H(:,:,tot)*kron(sig(:,:,ind),conj(sig(:,:,jnd))));
%         end
%     end
% end