function [Mout,strOut,Mcorr,alignmentDep] = alignToCentralBin(Min,mask,rcIn,visualize)
% Mout = alignToCentralBin(Min,mask)
% Min is a 9 x Nx x NAlines x Nbins SO(3) symmetric (linear) rotation 
% matrix. alignToCentralBin estimates the required correction to align all
% spectral bins to the central bin. Only points that are true within mask
% are used.
% Additional input arguments are rc(3 x Nbins), the correction vector, 
% impose a defined alignment. The local optimized rc is still computed for
% book-keeping.
% 'visualize' triggers the display of the correction.
% Optional output arguments are rc (3 x Nbins), the correction vector 
% applied to each spectral bin, and alignErr, the error between the origanl
% bins and the averaged corrected SO(3) matrices.
%
% Modified Jan 26 2018, to improve the angle estimation.

if nargin<3 || isempty(rcIn)
    rcIn = [];
end
if nargin<4
    visualize = false;
end

dim = size(Min);
if numel(dim)<=3% manage exception of no spectral binning
    dim(4) = 1;
end
indwc = ceil(dim(4)/2);

% number of points to consider for the optimization
% sample these randomly within the B-scans in regions with sufficient DOP
% (mask)
Npoints = dim(2);% use roughly one point per A-line
inds = find(mask);
inds = inds(round(rand(Npoints,1)*(numel(inds)-1))+1);
[zzind,xxind] = ind2sub([dim(2),dim(3)],inds(:)');
lInd = sub2ind(size(Min),repmat((1:9)',1,Npoints),repmat(zzind,9,1),repmat(xxind,9,1));
lInd = bsxfun(@plus,lInd,shiftdim((0:dim(4)-1)*prod(dim(1:3)),-1));

% subsampled MM
Mtemp = double(Min(lInd));

trRef = sum(Mtemp([1,5,9],:,indwc),1);% trace of the central bin
OAref = decomposeRot(Mtemp(:,:,indwc));% optic axis orientation of central bin
phiref = atan2(OAref(2,:),OAref(1,:));
%    [~,~,PAref] = decomposeRot(Mtemp(:,:,indwc));% optic axis orientation of central bin

% search for the best rotation, starting from the central bin, and using
% the previous result as initial condition
xopt = [0;0];
rr = zeros(3,dim(4));
R = zeros(9,dim(4));
options = optimoptions('fminunc','Display','off');
for indw = [indwc+1:dim(4),indwc-1:-1:1]
    fun = @(x) trRetError(x,trRef,Mtemp(:,:,indw));
    xopt = fminunc(fun,xopt,options);
    rr(:,indw) = [xopt(1);xopt(2);0];
    R = makeRot(rr(:,indw)/2);
    loc = MatrixMultiply(R,MatrixMultiply(Mtemp(:,:,indw),R));

    OA = decomposeRot(loc);
    phi = atan2(OA(2,:),OA(1,:));
    rotEst(indw) = median(mod(phiref-phi+pi/2,pi)-pi/2);
%        [~,~,PA] = decomposeRot(loc);
%        rotEst(indw) = mean(real(asin(squeeze(bsxfun(@times,PA(1,:),PAref(2,:)) - bsxfun(@times,PA(2,:),PAref(1,:))))));

    C = makeRot([0;0;-rotEst(indw)]);
    rc(:,indw) = decomposeRot(MatrixMultiply(R,C));

    if indw == dim(4)% reset initial point
        xopt = [0;0];
    end
end

if isempty(rcIn)
    rcEff = rc;
else
    rcEff = rcIn;
end

for indw = 1:dim(4)
    if indw == indwc
        Mcorr(:,:,:,indwc) = Min(:,:,:,indwc);
    else
        RC = makeRot(rcEff(:,indw));
        RCT = RC([1,4,7,2,5,8,3,6,9]).*[1;1;-1;1;1;-1;-1;-1;1];
        % more efficient implementation of pre- and post-multiplying with
        % the same matrix
        RR = kron(reshape(RC,[3,3]).',reshape(RCT,[3,3]));
        Mcorr(:,:,:,indw) = reshape(RR*reshape(Min(:,:,:,indw),[9,dim(2)*dim(3)]),[9,dim(2),dim(3)]);
%        Mcorr(:,:,:,indw) = MatrixMultiply(RCT,MatrixMultiply(Min(:,:,:,indw),RC));
    end
end

Mout = mean(Mcorr,4);
if dim(4)>1
%    [Mout,P] = euclideanRotation(Mout);
%    % compute 'depolarization' of this averaging P are the eigenvalues of
%    % the matirx Mout'*Mout, i.e. the transpose product. the depolarization
%    % index of the SO3 matrix is sqrt(trace(Mout'*Mout)/3) = sqrt(sum(P)/3)
%    alignmentDep = squeeze(sqrt(sum(P,1)/3));
    Mout = euclideanRotation(Mout);
else 
    alignmentDep = ones(dim(2:3));
end


if nargout>1
    % also compute alignErr
    % sum of (Mcorr - Mout).^2 wherever mask is true, divided by number of
    % true elements in mask
    strOut.alignErrInit = squeeze(sum(sum(sum(bsxfun(@times,((bsxfun(@minus,Min,Mout)).^2),shiftdim(mask,-1)),3),2),1))/sum(sum(mask));    
    strOut.alignErrEff = squeeze(sum(sum(sum(bsxfun(@times,((bsxfun(@minus,Mcorr,Mout)).^2),shiftdim(mask,-1)),3),2),1))/sum(sum(mask));    
    strOut.rc = rc;
    strOut.rcEff = rcEff;
end

if visualize
    % sanity check    
    [~,ret,PA] = decomposeRot(Min);
    [~,cret,cPA] = decomposeRot(Mcorr);
    [~,aret,aPA] = decomposeRot(Mout);

    angleOrig = real(asin(squeeze(bsxfun(@times,PA(1,:,:,:),PA(2,:,:,indwc)) - bsxfun(@times,PA(2,:,:,:),PA(1,:,:,indwc)))));
    angleCorr = real(asin(squeeze(bsxfun(@times,cPA(1,:,:,:),cPA(2,:,:,indwc)) - bsxfun(@times,cPA(2,:,:,:),cPA(1,:,:,indwc)))));

    figure(2)
    clf
    subplot(2,3,1)
    imagesc(squeeze(ret(1,:,:)))
    title('Orig cum ret')

    subplot(2,3,4)
    imagesc(squeeze(cret(1,:,:)))
    title('Corrected cum ret')

    subplot(2,3,2)
    imagesc(angleOrig(:,:))
    title('Angle to central bin')

    subplot(2,3,5)
    imagesc(angleCorr(:,:))
    title('Corrected angle to central bin')

    subplot(2,3,3)
    imagesc(squeeze(aret))
    title('Mean cum ret')

    subplot(2,3,6)
    hold on
    plot(rr')
    plot(rotEst)

end

function out = trRetError(x,trTarget,M)
R = makeRot([x(1);x(2);0]);
Mtemp = MatrixMultiply(M,R(:));
out = sum((trTarget(:).' - sum(Mtemp([1,5,9],:),1)).^2);