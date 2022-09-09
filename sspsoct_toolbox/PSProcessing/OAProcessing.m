function out = OAProcessing(S1,S2,procStruct,outputLevel)
% out = OAProcessing(S1,S2,procStruct,outputLevel)
% Depth-resolved optic axis processing, without compensation for dynamic
% system transmission. Reconstruction is achieved non-recursively. 
% S1, S2 are the spectrally binned Stokes vectors.

% ouputLevel defines the level of detail contained in out
% outputLevel = 0 : only final outputs
% outputLevel = 1 : additional computation of intermediate results
%
% Martin Villiger, March 27, 2020
%

if nargin<4 || isempty(outputLevel)
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
cumulative = false;
unwrapMask = [];

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
    elseif strcmp(fnames{ind},'cumulative')
        cumulative = procStruct.cumulative;
    elseif strcmp(fnames{ind},'unwrapMask')
        unwrapMask = procStruct.unwrapMask;
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

% final normalization
QUVf1 = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
QUVf2 = sqrt(dot(S2f(:,:,:,2:4),S2f(:,:,:,2:4),4));

dop1 = QUVf1./S1f(:,:,:,1);
dop2 = QUVf2./S2f(:,:,:,1);
dop = 1/2*mean(dop1,3) + 1/2*mean(dop2,3);
  
S1n = S1f(:,:,:,2:4)./QUVf1;
S2n = S2f(:,:,:,2:4)./QUVf2;

intf = mean(QUVf1+QUVf2,3);% filtered intensity

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

% Symmetrize matrices and align between specral bins
if outputLevel < 1
    [MM,strOut] = compensateSystem(MM,dopm,dopTh,wcorr,rc);
else
    [MM,strOut,MMcorr] = compensateSystem(MM,dopm,dopTh,wcorr,rc);
    out.errOpt = strOut.errOpt;
    out.errInit = strOut.errInit;
    out.errEff = strOut.errEff;
    out.alignErr = strOut.alignErrEff;
    out.alignErrInit = strOut.alignErrInit;
    if outputLevel > 1
        out.OAmean = decomposeRot(MM);
        out.OAaligned = decomposeRot(MMcorr);
    end
end
out.wcorr = strOut.wcorr;
out.wcorrEff = strOut.wcorrEff;
out.rc = strOut.rc;
out.rcEff = strOut.rcEff;

%%


[sqMMinv,sqWinv,unwrapMask] = unwrapOAx(MM,dop,dopTh,unwrapMask);
dmn = MatrixMultiply(sqMMinv(:,1:end-1,:),MatrixMultiply(MM(:,2:end,:),sqMMinv(:,1:end-1,:)));
locoa = decomposeRot(dmn);

% construct V-correction term
VV = MatrixMultiply(makeRot(locoa/2),MatrixMultiply(sqMMinv([1,4,7,2,5,8,3,6,9],1:end-1,:),sqMMinv(:,2:end,:)));
vv = decomposeRot(VV);

% first, use DOP to mask meaningless areas
mask = dop>dopTh;
mask = mask(1:end-1,:)&mask(2:end,:);
%    maskDop = mask;

%    mask = imerode(mask,ones(3*fwz,1));

% clean this mask up with an additional criterion on the uniformity of
% the local optic axis over an axial range.
maskoa = abs(imfilter(exp(1i*squeeze(atan2(locoa(1,:,:),locoa(2,:,:)))),ones(fwz,1)/fwz))>.9;

mask = mask&maskoa;

vvphi = squeeze(vv(3,:,:));
vvphi(~mask) = 0;    
vvcorr = shiftdim(cumsum(cat(1,zeros(1,size(vvphi,2)),vvphi(1:end-1,:))),-1);

oa(1,:,:) = locoa(1,:,:).*cos(vvcorr) - locoa(2,:,:).*sin(vvcorr);
oa(2,:,:) = locoa(1,:,:).*sin(vvcorr) + locoa(2,:,:).*cos(vvcorr);
oa(3,:,:) = zeros(size(locoa(1,:,:)));

locoa = cat(2,zeros(3,1,size(locoa,3)),locoa);
oa = cat(2,zeros(3,1,size(oa,3)),oa);    
% if outputLevel>0
%     % without depth correction
%     Omegaf = imfilter(permute(locoa,[2,3,1]),ones(fwz,1)/fwz);
%     retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
%     phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
%     out.phi2 = phi2;
%     out.ret2 = retFinal;
% end

if outputLevel>0
    if ~cumulative % without depth correction
        Omegaf = imfilter(permute(locoa,[2,3,1]),ones(fwz,1)/fwz);
        retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
        phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
        out.phi2 = phi2;
        out.ret2 = retFinal;
    else % outputting the cumulative signal
        W = decomposeRot(MM);
        Omegaf = permute(W,[2,3,1]);
        retFinal = sqrt(sum(Omegaf.^2,3))/pi*100;% scale pi to 100
        phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
        out.phi2 = phi2;
        out.ret2 = retFinal;
    end
end

Omegaf = imfilter(permute(oa,[2,3,1]),ones(fwz,1)/fwz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));


out.dop = dop;
out.mask = mask;
out.ret = retFinal;
out.phi = phi;

out.fwz = fwz;
out.fwx = fwx;
out.unwrapMask = unwrapMask;


