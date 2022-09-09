function [MM,outStruct,MMcorr] = compensateSystem(MM,dop,dopTh,wcorrIn,rcIn)
%[MMcorr] = compensateSystem(MM,dop,dopTh,wcorr,rcorr) returns the 
%symmetrized and system compensated SO3 rotation matrices, averaged across
%the spectral bins.
%Size of MM is [9,Nz,NAlines,Nbins], and each 3x3 matrix is a pure
%rotation matrix. dop has a size of [Nz,NAlines]. dopTh is a scalar 
%treshold for masking dop values when aligning the spectral bins.
%
%If wcorrIn (the symmetrization correction) and rcIn (bin alignment) are
%provided, they are direclty applied, rather than computing the specific
%correction vectors of the current data set. However, the symmetrization
%correction is still computed for book-keeping. 
%outStruct contains the 
%
%
%
% Martin Villiger March 2020


dim = size(MM);
if numel(dim)<=3 % manage exception of no spectral binning
    dim(4) = 1;
end

if nargin<4 || (~prod(size(wcorrIn) == [3,dim(4)])==1 && ~prod(size(rcIn) == [3,dim(4)])==1)
    wcorrIn = [];
    rcIn = [];
end


% Compute the correction necessary for making the system symmetric. Because
% this is quite fast, compute it anyway, even if wcorrIn is provided, for
% book-keeping over an entire volume. However, if wcorr and rcIn are
% provided, don't apply the correction, yet.

% Symmetrize matrices
if isempty(wcorrIn)
    [MM,strOut] = makeSymmetric(MM,dop);
else
    % don't apply correction, only compute errors
    [~,strOut] = makeSymmetric(MM,dop,wcorrIn,true);
end
outStruct.errOpt = strOut.errOpt;
outStruct.errInit = strOut.errInit;
outStruct.errEff = strOut.errEff;
outStruct.wcorr = strOut.wcorr;
outStruct.wcorrEff = strOut.wcorrEff;


% The bin-alignment is more computationally costly. Only perform it if rcIn
% is not available.
if isempty(rcIn)
%    if outputLevel < 1
    if nargin>2
        [MM,strOut,MMcorr] = alignToCentralBin(MM,dop>dopTh);
    else
        [MM,strOut] = alignToCentralBin(MM,dop>dopTh);
    end
    %    else
%        [MM,strOut,MMcorr] = alignToCentralBin(MM,dop>dopTh);
%        if outputLevel > 1
%            outStruct.OAmean = decomposeRot(MM);
%            outStruct.OAaligned = decomposeRot(MMcorr);
%        end
%    end
    outStruct.alignErrEff = strOut.alignErrEff;
    outStruct.alignErrInit = strOut.alignErrInit;
    outStruct.rc = strOut.rc;
    outStruct.rcEff = strOut.rcEff;
else
    % in this case, MM has not yet been symmetrized, and we now have to
    % apply to combined symmetry and alignment correction, and suppress
    % the v-components
    MMcorr = MM;
    for indw = 1:dim(4)
        Rsym = makeRot(wcorrIn(:,indw));
        Ralign = makeRot(rcIn(:,indw));
        RalignT = Ralign([1,4,7,2,5,8,3,6,9]).*[1;1;-1;1;1;-1;-1;-1;1];
        % What we need is RalignT*Rsym*MM*Ralign; we use A*B*C =
        % kron(C.',A)*vec(B)
        RR = kron(reshape(Ralign,[3,3]).',reshape(RalignT,[3,3])*reshape(Rsym,[3,3]));
        MMcorr(:,:,:,indw) = reshape(RR*reshape(MM(:,:,:,indw),[9,dim(2)*dim(3)]),[9,dim(2),dim(3)]);
    end
    
    MM = mean(MMcorr,4);
    MM = euclideanRotation(MM + MM([1,4,7,2,5,8,3,6,9]',:,:).*[1;1;-1;1;1;-1;-1;-1;1],true);
    
    mask = dop>dopTh;
    outStruct.alignErrInit = [];
    outStruct.alignErrEff = squeeze(sum(sum(sum(bsxfun(@times,((bsxfun(@minus,MMcorr,MM)).^2),shiftdim(mask,-1)),3),2),1))/sum(sum(mask));    
    outStruct.rc = [];
    outStruct.rcEff = rcIn;
end

