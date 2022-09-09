function [invM,detM] = matrix3inv(M)
    % computes the inverse of 3x3 matrices, assuming matrices linearized along first
    % dimension, in column order
    dim = size(M);
    
    if dim(1) ~= 9
        error('first dimnesion of the input should be 9');
    end
    M = reshape(M,9,[]);
    invM = M;
    invM(1,:) = M(5,:).*M(9,:)-M(8,:).*M(6,:);
    invM(2,:) = -(M(2,:).*M(9,:)-M(8,:).*M(3,:));
    invM(3,:) = (M(2,:).*M(6,:)-M(5,:).*M(3,:));
    invM(4,:) = -(M(4,:).*M(9,:)-M(7,:).*M(6,:));
    invM(5,:) = M(1,:).*M(9,:)-M(7,:).*M(3,:);
    invM(6,:) = -(M(1,:).*M(6,:)-M(4,:).*M(3,:));
    invM(7,:) = M(4,:).*M(8,:)-M(7,:).*M(5,:);
    invM(8,:) = -(M(1,:).*M(8,:)-M(7,:).*M(2,:));
    invM(9,:) = M(1,:).*M(5,:)-M(4,:).*M(2,:);
    
    detM = M(1,:).*invM(1,:)+M(4,:).*invM(2,:)+M(7,:).*invM(3,:);
    detMiv = detM;
    %detMiv(detM<1e-6) = 1;
    invM = invM./detMiv;
    invM = reshape(invM,dim);
    dim(1)=1;
    detM = reshape(detM,dim);
end

