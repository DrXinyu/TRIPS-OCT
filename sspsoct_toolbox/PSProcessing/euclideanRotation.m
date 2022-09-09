function [Rout,P,VV,err1,err2] = euclideanRotation(Min,boolDSym)
%euclideanRotationFast(Min) decomposes the 3x3 rotation matrix Min, in linear
%form (9 consecutive elements along 1st dimension) into the corresponding
%pure rotation matrix Rout, by separating it from a positive semidefinite
%component. Rather than using an iterative eigendecomposition of each
%matrix, we use Moakher, M. Means and Averaging in the Group of Rotations. 
% SIAM Journal on Matrix Analysis and Applications (2002). In short, the
% arithmetic mean of all rotation matrices consists of the mean rotation
% matrix we are seeking, and an additional, depolarization matrix, which is
% positive-semidefinite: Rm = R*V. To isolate R, we compute R =
% Rm*inv(sqrtm(Rm'*Rm)). To efficiently compute the eigenvalues and
% eigenvectors of Rm'*Rm, the formula by Smith, O. K. Eigenvalues of a 
% symmetric 3 × 3 matrix. Commun. ACM 4, 168 (1961), is employed.
% Also see https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
% 
% When averaging D-symmetric sums of a matrix, i.e. M + D*M.'*D, there
% arises the possibility that when taking the matrix root, rather than
% taking the positive square roots of the eigenvalues, we should take the
% negative of the two smaller eigenvalues. This case is managed when
% boolDSym is set to true. Default is false. 
% 
% Optional output arguments are P, the computed eigenvalues of the
% positive-semitdefinite matrix, its eigenvectors VV, and err1 and err2,
% the squared 'distance' from [-1,-1,1] of the final rotation matrix.

% make sure input is double
tt = whos('Min');
Min = double(Min);

if nargin<2
    boolDSym = false;
end

% Compute Rm'*Rm
dim = [size(Min),1];
R2 = MatrixMultiply(Min([1,4,7,2,5,8,3,6,9],:),Min(:,:));

% conventional determininant of Min
ddet = Min(1,:).*Min(5,:).*Min(9,:) + Min(2,:).*Min(6,:).*Min(7,:) + Min(3,:).*Min(4,:).*Min(8,:) - Min(3,:).*Min(5,:).*Min(7,:) - Min(1,:).*Min(6,:).*Min(8,:) - Min(2,:).*Min(4,:).*Min(9,:);



% computation of the eigenvalues of R2 as outlined in https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
p1 = R2(4,:).^2 + R2(7,:).^2 + R2(8,:).^2;

q = sum(R2([1,5,9],:),1)/3;
p2 = (R2(1,:) - q).^2 + (R2(5,:) - q).^2 + (R2(9,:) - q).^2 + 2 * p1;
p = sqrt(p2 / 6);
R2p = R2;
for ind = [1,5,9]
    R2p(ind,:) = R2(ind,:) - q;
end
B = bsxfun(@rdivide,R2p,p);
r = (B(1,:).*B(5,:).*B(9,:) + B(2,:).*B(6,:).*B(7,:) + B(3,:).*B(4,:).*B(8,:) - B(3,:).*B(5,:).*B(7,:) - B(1,:).*B(6,:).*B(8,:) - B(2,:).*B(4,:).*B(9,:))/2;
r = min(max(r,-1),1);
phi = acos(r) / 3;

% the eigenvalues satisfy eig1 <= eig2 <= eig3
P(3,:) = q + 2 .* p .* cos(phi);
P(1,:) = q + 2 .* p .* cos(phi + (2*pi/3));
P(2,:) = 3 .* q - P(1,:) - P(3,:);     % since trace(A) = eig1 + eig2 + eig3

P = abs(P);% ensure all eigenvalues are positive

% computation of the eigenvectors of R2
AA = cat(1,sum(bsxfun(@times,R2([1,4,7],:),R2(1:3,:)),1),sum(bsxfun(@times,R2([2,5,8],:),R2(1:3,:)),1),sum(bsxfun(@times,R2([3,6,9],:),R2(1:3,:)),1));
v1 = AA - bsxfun(@times,R2(1:3,:),P(2,:)+P(3,:));
v1(1,:) = v1(1,:) + P(2,:).*P(3,:);
v2 = AA - bsxfun(@times,R2(1:3,:),P(1,:)+P(3,:));
v2(1,:) = v2(1,:) + P(1,:).*P(3,:);
v3 = AA - bsxfun(@times,R2(1:3,:),P(1,:)+P(2,:));
v3(1,:) = v3(1,:) + P(1,:).*P(2,:);

%If the input matrix already is a pure retarder the three eigenvalues are
%1, and v1, v2, or v3 may be 0. In this case, v1 and v3 can be initalized to
%a random value. v2 has to be different from either.
mask = sum(v1==0,1)==3;
if sum(mask)>0
    v1(:,mask) = randn(3,gather(sum(mask)));
end
mask = sum(v3==0,1)==3;
if sum(mask)>0
    v3(:,mask) = randn(3,gather(sum(mask)));
end
% re-initialize v2 in all cases
v2 = cross(v1,v3,1);
% if v1 and v3 are collinear, this results in 0
mask = sum(v2==0,1)==3;
if sum(mask)>0
    [~,inds] = min(abs(v1(:,mask)),[],1);
    zz = zeros(3,numel(inds));
    zz(inds + (0:numel(inds)-1)*3) = 1;
    v2(:,mask) = zz;
end

V1 = bsxfun(@rdivide,v1,sqrt(sum(v1.^2,1)));
V2 = bsxfun(@rdivide,v2,sqrt(sum(v2.^2,1)));
V3 = bsxfun(@rdivide,v3,sqrt(sum(v3.^2,1)));

% find which eigenvalues are further apart
lInds = sign((P(2,:)-P(1,:))-(P(3,:)-P(2,:)));

% Force the eigenvectors to be orthogonal (which they should be, but maybe
% not precisely are). 
% The most reliable eigenvector is that one originating of the two more
% closely spaced eigenvalues;
% lInds<=0 ;i.e. 1 and 2 are close, V3 is reliable
V1c = V1;
V2c = V2;
V3c = V3;
% correct V2c and V1c
V2c(:,lInds<=0) = V2(:,lInds<=0) - bsxfun(@times,V3(:,lInds<=0),sum(V3(:,lInds<=0).*V2(:,lInds<=0),1));
V1c(:,lInds<=0) = cross(V3(:,lInds<=0),V2c(:,lInds<=0),1);

% correct V2c and V3c
V2c(:,lInds>0) = V2(:,lInds>0) - bsxfun(@times,V1(:,lInds>0),sum(V1(:,lInds>0).*V2(:,lInds>0),1));
V3c(:,lInds>0) = cross(V1(:,lInds>0),V2c(:,lInds>0),1);
% check that Vn's are orthogonal

V1c = bsxfun(@rdivide,V1c,sqrt(sum(V1c.^2,1)));
V2c = bsxfun(@rdivide,V2c,sqrt(sum(V2c.^2,1)));
V3c = bsxfun(@rdivide,V3c,sqrt(sum(V3c.^2,1)));
 
VV = cat(1,V1c,V2c,V3c);

VV(:,p1==0) = repmat(reshape(eye(3),[9,1]),1,sum(p1==0));

% compute Rm/sqrtm(Rm'*Rm)
Rout = MatrixMultiply(Min(:,:),VV);
Rout(1:3,:) = bsxfun(@rdivide,Rout(1:3,:),sqrt(P(1,:)).*sign(ddet)); % sign(ddet) ensured that det(Rout)=+1
Rout(4:6,:) = bsxfun(@rdivide,Rout(4:6,:),sqrt(P(2,:)));
Rout(7:9,:) = bsxfun(@rdivide,Rout(7:9,:),sqrt(P(3,:)));

Rout = MatrixMultiply(Rout,VV([1,4,7,2,5,8,3,6,9],:));

if boolDSym% manage trivial solution diag([-1,-1,1])
    RRout = MatrixMultiply(Min(:,:),VV);
    RRout(1:3,:) = bsxfun(@rdivide,RRout(1:3,:),-sqrt(P(1,:)).*sign(ddet)); % sign(ddet) ensured that det(Rout)=+1
    RRout(4:6,:) = bsxfun(@rdivide,RRout(4:6,:),-sqrt(P(2,:)));
    RRout(7:9,:) = bsxfun(@rdivide,RRout(7:9,:),sqrt(P(3,:)));
    RRout = MatrixMultiply(RRout,VV([1,4,7,2,5,8,3,6,9],:));
    
    % compare which diagonal is further from [-1,-1,1]
    err1 = reshape(sum((Rout([1,5,9],:)-[-1;-1;1]).^2,1),dim(2:end));
    err2 = reshape(sum((RRout([1,5,9],:)-[-1;-1;1]).^2,1),dim(2:end));
    mask = 10*sum((Rout([1,5,9],:)-[-1;-1;1]).^2,1) < sum((RRout([1,5,9],:)-[-1;-1;1]).^2,1);% factor of 10 was found empirically
    Rout(:,mask) = RRout(:,mask);
else
    err1 = [];
    err2 = [];
end

Rout = reshape(Rout,dim);
%Rout = cast(Rout,tt.class);

if nargout>1
    P = reshape(P,[3,dim(2:end)]);
end
if nargout>2
    VV = reshape(VV,dim);
end