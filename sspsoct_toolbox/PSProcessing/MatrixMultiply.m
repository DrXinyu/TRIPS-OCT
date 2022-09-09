function  Out = MatrixMultiply(A,B)
% computes matrix products, assuming matrices linearized along first
% dimension, in column order

N = sqrt(size(A,1));% number of elements along one side of the matrix

dimA = size(A);
dimB = size(B);

dimA = [dimA,ones(1,numel(dimB)-numel(dimA))];
dimB = [dimB,ones(1,numel(dimA)-numel(dimB))];

S1.type = '()';
S1.subs = repmat({':'},1,ndims(A));
S2.type = '()';
S2.subs = repmat({':'},1,ndims(B));
 
Out = zeros(max(dimA,dimB));
for ind3 = 1:N
    ind1 = repmat((ind3-1)*N+1:ind3*N,[1,N])';
    ind2 = repmat((ind3:N:N*N),[N,1]);
    S1.subs{1} = ind1;
    S2.subs{1} = ind2;
    Out = Out + subsref(A,S1).*subsref(B,S2);
end

% ind3 = 1;
% ind1 = repmat((ind3-1)*N+1:ind3*N,[1,N])';
% ind2 = repmat((ind3:N:N*N),[N,1]);
% Out(:,:) = bsxfun(@times,A(ind1,:),B(ind2,:));
% for ind3 = 2:N
%     ind1 = repmat((ind3-1)*N+1:ind3*N,[1,N])';
%     ind2 = repmat((ind3:N:N*N),[N,1]);
%     Out(:,:) = Out(:,:) + bsxfun(@times,A(ind1,:),B(ind2(:),:));
% end
% 
% Out = reshape(Out,dim);
