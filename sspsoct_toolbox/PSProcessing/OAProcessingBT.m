function out = OAProcessingBT(S1,S2,int,procStruct,outputLevel)
% out = OAProcessingBT(S1,S2,int,procStruct,outputLevel)
% Depth-resolved optix axis in bench-top configuration
% S1, S2 are the spectrally binned Stokes vectors.

% ouputLevel defines the level of detail contained in out
% outputLevel = 0 : only final outputs
% outputLevel = 1 : additional computation of intermediate results
%
% Martin Villiger, June 21, 2017
%

if nargin<5 || isempty(outputLevel)
    outputLevel = 0;% default
end

% parameters
roiz = 1:1000; % axial range to consider
fwx = 12;% lateral filtering
fwaxial = 1;% for comparison without spectral binning
dopTh = .7; % for suppressing meaningless regions
fwz = 5; % axial filtering range of local ret
dzres = 4.8;
rc = [];
wcorr = [];
boolFastSym = true; % fast computation in makeSymmetric, setting V-component to zero

fnames = fieldnames(procStruct);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'fwx')
        fwx = procStruct.fwx;
    elseif strcmp(fnames{ind},'fwz')
        fwz = procStruct.fwz;
    elseif strcmp(fnames{ind},'roiz')
        roiz = procStruct.roiz;
    elseif strcmp(fnames{ind},'rc')
        rc = procStruct.rc;
    elseif strcmp(fnames{ind},'wcorr')
        wcorr = procStruct.wcorr;
    elseif strcmp(fnames{ind},'dzres')
        dzres = procStruct.dzres;
    elseif strcmp(fnames{ind},'dopTh')
        dopTh = procStruct.dopTh;
    elseif strcmp(fnames{ind},'fwaxial')
        fwaxial = procStruct.fwaxial;
    elseif strcmp(fnames{ind},'boolFastSym')
        boolFastSym = procStruct.boolFastSym;
    end
end

dim = size(S1);
% manage case of no spectral binning
if numel(dim)<4 % no spectral binning used
    S1 = permute(S1,[1,2,4,3]);% introduce third dimension
    S2 = permute(S2,[1,2,4,3]);% introduce third dimension
    dim = size(S1);
end

if dim(4) == 3 % 3-component Stokes vector was provided
%    I1 = sqrt(dot(S1,S1,4));
%    I2 = sqrt(dot(S2,S2,4));
    S1 = cat(4,sqrt(dot(S1,S1,4)),S1);
    S2 = cat(4,sqrt(dot(S2,S2,4)),S2);
end

nx = (round(fwx*1.5)-1)/2;
nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
h = exp(-nx.^2);
if fwaxial>1
    nz = (round(fwaxial*1.5)-1)/2;
    nz = linspace(-nz,nz,round(fwaxial*1.5))*2*sqrt(log(2))/fwaxial;
    h = exp(-nz(:).^2)*h;
end
h = h/sum(h(:));
S1f = imfilter(S1,h,'circular');
S2f = imfilter(S2,h,'circular');
%I1f = imfilter(I1,h,'circular');
%I2f = imfilter(I2,h,'circular');

% final normalization
QUVf1 = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
QUVf2 = sqrt(dot(S2f(:,:,:,2:4),S2f(:,:,:,2:4),4));

dop1 = QUVf1./S1f(:,:,:,1);
dop2 = QUVf2./S2f(:,:,:,1);
% dop1 = QUVf1./I1f;
% dop2 = QUVf2./I2f;
dop = 1/2*mean(dop1,3) + 1/2*mean(dop2,3);
  
S1n = S1f(:,:,:,2:4)./QUVf1;
S2n = S2f(:,:,:,2:4)./QUVf2;

% construct orthonormal tripod for these data points
na = S1n + S2n;
nb = S1n - S2n;
na = na./repmat(sqrt(dot(na,na,4)),[1,1,1,3]);
nb = nb./repmat(sqrt(dot(nb,nb,4)),[1,1,1,3]);
S1n = (na + nb)/sqrt(2);
S2n = (na - nb)/sqrt(2);

% construct 3x3 rotation matrix at each pixel
S3n = cross(S1n,S2n,4);

% generate measurement matrix; dimensions: 3x3,Nz,NAlines,Nw
MM = permute(cat(4,S1n,S2n,S3n),[4,1,2,3]);

if outputLevel > 1
    out.OA1 = decomposeRot(MM);
end

%%
dopm = dop;
dopm(roiz(end):end,:) = 0;

% Symmetrize matrices
if outputLevel < 1
    [MM,strOut] = makeSymmetric(MM,dopm,wcorr,[],boolFastSym);
else
    [MM,strOut] = makeSymmetric(MM,dopm,wcorr,[],boolFastSym);
    out.errOpt = strOut.errOpt;
    out.errInit = strOut.errInit;
    out.errEff = strOut.errEff;
end
out.wcorr = strOut.wcorr;
out.wcorrEff = strOut.wcorrEff;
%%
[MM,strOut] = alignToCentralBin(MM,dop>dopTh,rc);
out.alignErr = strOut.alignErrEff;
out.alignErrInit = strOut.alignErrInit;
out.rc = strOut.rc;
out.rcEff = strOut.rcEff;
    
%%
% generate mask to exclude points with low DOP from correcting the local
% retardation
mmask = (dop>dopTh);
mmask = imdilate(imerode(mmask,ones(fwz,1)),ones(fwz,1));

dim = size(MM);
w = zeros([3,dim(2),dim(3)]);
if outputLevel>0
    w2 = zeros([3,dim(2),dim(3)]);
end

N = repmat([1;0;0;0;1;0;0;0;1],[1,1,dim(3)]);
%N = permute(makeRot(cat(1,zeros(2,numel(Vcorr)),-Vcorr)),[1,3,2]);

if outputLevel>0
    W = decomposeRot(MM);
    Msqinv = makeRot(-W/2);
end
for indz = 2:roiz(end)%refInds(3):roiz(end)
    
    % D*N(indz)*D*Mtot(indz+1)*N'(indz)
    nloc = MatrixMultiply(bsxfun(@times,N,[1;1;-1;1;1;-1;-1;-1;1]),MatrixMultiply(MM(:,indz,:,:),N([1,4,7,2,5,8,3,6,9],:,:,:)));

    Wloc = decomposeRot(nloc);
 
    w(:,indz,:) = squeeze(Wloc);
    nlocroot = makeRot(Wloc/2);

    Nnext = MatrixMultiply(nlocroot,N);
    N(:,:,mmask(indz,:)>0) = Nnext(:,:,mmask(indz,:)>0);

    if outputLevel>0
        % without depth correction
        nloc = MatrixMultiply(Msqinv(:,indz-1,:,:),MatrixMultiply(MM(:,indz,:,:),Msqinv(:,indz-1,:,:)));
%        nloc = MatrixMultiply(MM(:,indz,:,:),MM([1,4,7,2,5,8,3,6,9],indz-1,:,:));
        Wloc = decomposeRot(nloc);
        w2(:,indz,:) = squeeze(Wloc);
    end
end

if outputLevel>0
    % without depth correction
    Omegaf = imfilter(permute(w2,[2,3,1]),ones(fwz,1)/fwz);
    retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
    phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
    out.phi2 = phi2;
    out.ret2 = retFinal;
end
    
Omegaf = imfilter(permute(w,[2,3,1]),ones(fwz,1)/fwz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
          
out.dop = dop;
out.mask = mmask;
out.ret = retFinal;
out.phi = mod(phi + pi,2*pi)-pi;

out.fwz = fwz;
out.fwx = fwx;
out.int = int;


