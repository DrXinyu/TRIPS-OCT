function H = M2H(M)
% converts and input Mueller matrix M into its corresponding 'H'-matrix,
% following Cloude.

N = size(M,3);

A = 1/sqrt(2)*[1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0];

for tot = 1:N
    H = A'*M(:,:,tot)*A;
    H([2,4,6,8,9,11,13,15]) = H([9,11,13,15,2,4,6,8]);
end
    