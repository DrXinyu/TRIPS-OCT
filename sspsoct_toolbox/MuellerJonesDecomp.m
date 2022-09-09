function [r,d] = MuellerJonesDecomp(Min,polar)
% Performs simultaneous decomposition of a deterministic, i.e.
% Mueller-Jones matrix, which has only retardance and diattenuation.
% 
% The function converts the Mueller Jones matrix to its H-matrix. If M is
% indeed deterministic, then the H-matrix if of rank 1 and can be written
% as H = j*j', where j is the vectorize version of the underlying Jones
% matrix.

if nargin<2 
    polar = false;
end

dim = size(Min);
if dim(1) == 4 && dim(2) == 4
    if numel(dim)>2
        dim = cat(2,16,dim(3:end));
        Min = reshape(Min,cat(2,16,dim(2:end)));
    else
        dim = [16,1];
        Min = Min(:);
    end
end

% convert M to H matrix
A = 1/sqrt(2)*[1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0];
At = A';

H = MatrixMultiply(At(:),MatrixMultiply(Min,A(:)));
H([2,4,6,8,9,11,13,15],:) = H([9,11,13,15,2,4,6,8],:);
J = H([1;3;2;4],:)./sqrt(H(1,:));

[r,d] = JonesDecomp(J,polar);

r = reshape(r,[3,dim(2:end)]);
d = reshape(d,[3,dim(2:end)]);
