function out = PSProcessGlobalSymmetric(S1,S2,procStruct)
%out = PSProcessGlobalSymmetric(S1,S2,procStruct) - computes local 
%birefringence of the Stokes vectors S1 and S2, corresponding to the Stokes 
%vectors measured for an input polarization state modulated between states 
%orthogonal on the Poincaree sphere, using spectral binning. In contrast to 
%PSPRocess, which estimates the alignment between the spectral bins for 
%each A-line, PSProcessGlobal computes a constant alignment across the 
%entire B-scan. Furthermore, it aligns the spectral bins before computing
%the local, depth-resolved retardation.

% the only two mandatory arguments
fwx = procStruct.fwx;
dz = procStruct.dz;
fwz = 1;
dopTh = 0.6;
wcorr = [];
rc = [];
dzres = 4.8;

fnames = fieldnames(procStruct);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'fwz')
        fwz = procStruct.fwz;
    elseif strcmp(fnames{ind},'rc')
        rc = procStruct.rc;
    elseif strcmp(fnames{ind},'wcorr')
        wcorr = procStruct.wcorr;
    elseif strcmp(fnames{ind},'dzres')
        dzres = procStruct.dzres;
    elseif strcmp(fnames{ind},'dopTh')
        dopTh = procStruct.dopTh;
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
    S1 = cat(4,sqrt(dot(S1,S1,4)),S1);
    S2 = cat(4,sqrt(dot(S2,S2,4)),S2);
end

nx = (round(fwx*1.5)-1)/2;
nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
h = exp(-nx.^2);
if fwz>1
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

%%
dopm = dop;

% Symmetrize matrices and align between specral bins
[MM,strOut] = compensateSystem(MM,dopm,dopTh,wcorr,rc);
out.errOpt = strOut.errOpt;
out.errInit = strOut.errInit;
out.errEff = strOut.errEff;
out.alignErr = strOut.alignErrEff;
out.alignErrInit = strOut.alignErrInit;
out.wcorr = strOut.wcorr;
out.wcorrEff = strOut.wcorrEff;
out.rc = strOut.rc;
out.rcEff = strOut.rcEff;


dmn = MatrixMultiply(MM(:,2:end,:),MM([1,4,7,2,5,8,3,6,9],1:end-1,:));
locoa = cat(2,zeros(3,1,size(MM,3)),decomposeRot(dmn));

Omegaf = imfilter(permute(locoa,[2,3,1]),ones(dz,1)/dz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;

out.dop = dop;
out.ret = retFinal;
out.S1f = S1f;
out.S2f = S2f;



