function [MMcorr,outStruct,MMcorrraw] = makeSymmetric(MM,dop,wcorrIn,boolNoMMcorr,boolFastSym)
%[MMcorr] = makeSymmetric(MM,dop,wcorr) returns the symmetrized measurement 
%matrix. Size of MM is [9,Nz,NAlines,Nbins], and each 3x3 matrix is a pure
%rotation matrix. dop has a size of [Nz,NAlines]. The optional input
%argument wcorrIn results in directly applying the correction wcorrIn. The
%local wcorr is still computed, together with its specific error and
%returned as outStruct.wcorr.
% 
%MMcorr is the symmetrized MM, with the remaining component along the
%V-coordinate set to zero. If boolNoMMCorr is set to ture, its computation
%is suppressed. If boolFastSym is set to true (default), the V-component is
%indeed simply forced to zeros. if boolFastSym is false, the least square
%solution is copmuted using euclideanRotation(M + D*M.'*D), which is
%slower, but more accurate (in the lsq sense).
%
%outStruct contains fields wcorr, the computed wcorr, wcorrEff, the
%effectively used one (wcorrIn, if provided, otherwise wcorr). errInit,
%errOpt and errEff, the normalized V-coordinate components for each
%spectral bin, for the initial MM, the optimized symmetrization and the
%effectively computed one (different from opt if wcorrIn provided).
% MMcorrraw is the symmetrized MM, but without setting the V-coordinate to 
% zero.
%
% This code operates directly on the supplied Mueller matrices, without
% constructing the corresponding Jones matrices
%
% Martin Villiger November 26, 2017
% Modified March 2020: Apply wcorrIn, if provided but still compute the
% local wcorr.


dim = size(MM);
if numel(dim)<=3 % manage exception of no spectral binning
    dim(4) = 1;
end

if nargin<3 || ~prod(size(wcorrIn) == [3,dim(4)])==1
    wcorrIn = [];
end

if nargin<4 || isempty(boolNoMMcorr)
    boolNoMMcorr = false;
end

if nargin<5
    boolFastSym = true;
end

MMcorr = MM;
dopTh = 0.7;
ww = dop>dopTh;

% First, average the SO3 matrices over the mask ww
for ind = 1:dim(1)
    %size(MMproj) is 9 by number of bins
    MMproj(ind,:) = ww(:)'*reshape(MM(ind,:,:,:),numel(ww),dim(4));
end

A = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0]/sqrt(2);
for indw = 1:dim(4)
    % find estimation matrix, convert MMproj to H-matrix, then multiply
    % with diag([-1,1,-1,1]) on both sides
    F = conj(A'*[sum(ww(:)),0,0,0;cat(2,zeros(3,1),reshape(MMproj(:,indw),[3,3]))]*A);
    H = reshape(F([1,9,3,11,5,13,7,15,2,10,4,12,6,14,8,16]'),[4,4]); 
    [a,b] = eig(diag([-1,1,-1,1])*H*diag([-1,1,-1,1]));
    [~,pos] = min(abs(diag(b)));
    cc = a(:,pos);

    Jcorr = [cc(2),cc(4);cc(1),cc(3)];
    Jcorr = Jcorr/sqrt(det(Jcorr));

    [rr,dd] = JonesDecomp(Jcorr);% convert back to Stokes, ignoring diattanuation component

    % collect correction rotation vector for all spectral bins
    wcorr(:,indw) = real(rr);
    ddcorr(:,indw) = real(dd);

    % error computation
    jvec = [Jcorr(2,1),Jcorr(1,1),Jcorr(2,2),Jcorr(1,2)].';
    errOpt(indw) = real(jvec'*diag([-1,1,-1,1])*H*diag([-1,1,-1,1])*jvec/4/sum(ww(:))); % factor of 4 is as reference for entire energy; each pixel within ww is a normalized Jones matrix with sum(abs(J(:)).^2) = 2
    % more pragmatically, a pure V-component matrix M with pi
    % retardation results in an error of 1
    jvec = [0, 1, 1, 0]';
    errInit(indw) = real(jvec'*diag([-1,1,-1,1])*H*diag([-1,1,-1,1])*jvec/4/sum(ww(:)));


    if ~isempty(wcorrIn)
        Jcorr = makeJones(wcorrIn(:,indw));
        jvec = [Jcorr(2,1),Jcorr(1,1),Jcorr(2,2),Jcorr(1,2)].';
        errEff(indw) = real(jvec'*diag([-1,1,-1,1])*H*diag([-1,1,-1,1])*jvec/4/sum(ww(:)));
    end
    
end

% unwrap correction across spectral bins
ss = cat(2,0,mod(cumsum(sqrt(sum(diff(wcorr,[],2).^2,1))>pi),2));
if sum(ss)>dim(4)/2
    ss = 1-ss;
end
retcorr = sqrt(sum(wcorr.^2,1));
wcorr = wcorr - bsxfun(@times,bsxfun(@rdivide,wcorr,retcorr),ss*2*pi);

% maintain a retardation <pi
if mean(sqrt(sum(wcorr.^2,1)))>pi
    wcorr = wcorr - bsxfun(@rdivide,wcorr,sqrt(sum(wcorr.^2,1)))*2*pi;
end

% use the computed wcorr only of no wcorrIn has been provided
if isempty(wcorrIn)
    wcorrEff = wcorr;
    errEff = errOpt;
else
    wcorrEff = wcorrIn;
end

% Apply correction; only do this, if MM is requested as output argument
if ~boolNoMMcorr
    % attempt to speed this up
    for indw = 1:dim(4)
        Rcorr = makeRot(wcorrEff(:,indw));
        MMcorr(1,:,:,indw) = MM(1,:,:,indw)*Rcorr(1) +  MM(2,:,:,indw)*Rcorr(4) + MM(3,:,:,indw)*Rcorr(7);
        MMcorr(2,:,:,indw) = MM(1,:,:,indw)*Rcorr(2) +  MM(2,:,:,indw)*Rcorr(5) + MM(3,:,:,indw)*Rcorr(8);
        MMcorr(3,:,:,indw) = MM(1,:,:,indw)*Rcorr(3) +  MM(2,:,:,indw)*Rcorr(6) + MM(3,:,:,indw)*Rcorr(9);

        MMcorr(4,:,:,indw) = MM(4,:,:,indw)*Rcorr(1) +  MM(5,:,:,indw)*Rcorr(4) + MM(6,:,:,indw)*Rcorr(7);
        MMcorr(5,:,:,indw) = MM(4,:,:,indw)*Rcorr(2) +  MM(5,:,:,indw)*Rcorr(5) + MM(6,:,:,indw)*Rcorr(8);
        MMcorr(6,:,:,indw) = MM(4,:,:,indw)*Rcorr(3) +  MM(5,:,:,indw)*Rcorr(6) + MM(6,:,:,indw)*Rcorr(9);

        MMcorr(7,:,:,indw) = MM(7,:,:,indw)*Rcorr(1) +  MM(8,:,:,indw)*Rcorr(4) + MM(9,:,:,indw)*Rcorr(7);
        MMcorr(8,:,:,indw) = MM(7,:,:,indw)*Rcorr(2) +  MM(8,:,:,indw)*Rcorr(5) + MM(9,:,:,indw)*Rcorr(8);
        MMcorr(9,:,:,indw) = MM(7,:,:,indw)*Rcorr(3) +  MM(8,:,:,indw)*Rcorr(6) + MM(9,:,:,indw)*Rcorr(9);
    end

    if nargout>2
        % maintain MMcorr in memory before setting V-components to zero 
        MMcorrraw = MMcorr;
    end

    % set V-component to zero
    if boolFastSym
        OA = decomposeRot(MMcorr);
        OA(3,:) = 0;
        MMcorr = makeRot(OA);
    else
        MMcorr = euclideanRotation(MMcorr + MMcorr([1,4,7,2,5,8,3,6,9]',:,:,:).*[1;1;-1;1;1;-1;-1;-1;1],true);
    end
else
    % if boolNoMMcorr is true, simply assign empty matrix
    MMcorr = [];
    MMcorrraw = [];
end

outStruct.wcorr = wcorr;
outStruct.wcorrEff = wcorrEff;

outStruct.errOpt = errOpt;
outStruct.errInit = errInit;
outStruct.errEff = errEff;

