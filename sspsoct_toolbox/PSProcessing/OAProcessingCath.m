function out = OAProcessingCath(S1,S2,int,procStruct,outputLevel)
% out = OAProcessingCath(S1,S2,int,procStruct,outputLevel)
% Depth-resolved optix axis in catheter configuration
% S1, S2 are the spectrally binned Stokes vectors.

% ouputLevel defines the level of detail contained in out
% outputLevel = 0 : only final outputs
% outputLevel = 1 : additional info on processing steps
% outputLevel = 2 : intermediate optic axes results
%
% Martin Villiger, June 22, 2017
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
fw = 30; % lateral filtering to estimate surface signal and rotation due to sheath birefringence
dzres = 4.8;
rc = [];
wcorr = [];
cumulative = false; % flag for processing cumulative signal

fnames = fieldnames(procStruct);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'fwx')
        fwx = procStruct.fwx;
    elseif strcmp(fnames{ind},'fw')
        fw = procStruct.fw;
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
dop = 1/2*mean(dop1,3) + 1/2*mean(dop2,3);
  
S1n = S1f(:,:,:,2:4)./repmat(QUVf1,[1,1,1,3]);
S2n = S2f(:,:,:,2:4)./repmat(QUVf2,[1,1,1,3]);

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
    [MM,strOut] = makeSymmetric(MM,dopm,wcorr);
else
    [MM,strOut,MMraw] = makeSymmetric(MM,dopm,wcorr);
    out.errInit = strOut.errInit;
    out.errEff = strOut.errEff;
    if outputLevel > 1
        out.OA2raw = decomposeRot(MMraw);
        out.OA2 = decomposeRot(MM);
    end
end
out.errOpt = strOut.errOpt;
out.wcorr = strOut.wcorr;
out.wcorrEff = strOut.wcorrEff;

%%

if outputLevel < 1
    [MM,strOut] = alignToCentralBin(MM,dop>dopTh,rc);
else
    [MM,strOut,MMcorr] = alignToCentralBin(MM,dop>dopTh,rc);
    out.alignErr = strOut.alignErrEff;
    out.alignErrInit = strOut.alignErrInit;
    if outputLevel > 1
        out.OAmean = decomposeRot(MM);
        out.OAaligned = decomposeRot(MMcorr);
    end
end
out.rc = strOut.rc;
out.rcEff = strOut.rcEff;
    
%%
% detect catheter
cath = findCatheter(int);

if outputLevel > 0
    out.cath = cath;
end

% identify signal from catheter
%Itotf = I1f + I2f;

% Take intensity weighted average between first and second peak
mm = (bsxfun(@minus,(1:size(dop,1))',cath(1,:)-10)>0)&(bsxfun(@minus,(1:size(dop,1))',cath(2,:)+10)<0);
rrz = max(min(cath(1,:))-10,1):max(cath(2,:))+10;
mm = mm(rrz,:);
ww = shiftdim(repmat(max(dop(rrz,:)*2-1,0).*mm,[1,1,size(MM,4)]),-1);
Mm = squeeze(bsxfun(@rdivide,sum(bsxfun(@times,MM(:,rrz,:,:),ww),2),sum(ww*9,2)));
Mm = imfilter(Mm,ones(1,fw)/fw,'circular');
Mm = euclideanRotation(Mm);
WFit = decomposeRot(Mm);

if outputLevel > 0
    out.surf = WFit;
end

%% unwrap fit along the B-scan
ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(WFit,[],2).^2,1))>pi,2),2)),[3,1])>0;
delta = 2*pi*bsxfun(@rdivide,WFit,sqrt(sum(WFit.^2,1)));

WFitc = WFit;
WFitc(ss) = WFit(ss) - delta(ss);

% Check if overall a subtraction/addition of pi would make more sense
n = -1:1;
for ind = 1:3
    WW = WFitc + 2*pi*n(ind)*bsxfun(@rdivide,WFitc,sqrt(sum(WFitc.^2,1)));
    meanRet(ind) = mean(mean(sqrt(sum(WW.^2,1))));
end

[~,mp] = min(meanRet);
WW = WFitc + 2*pi*n(mp)*bsxfun(@rdivide,WFitc,sqrt(sum(WFitc.^2,1)));
WW(3,:,:) = 0;

if outputLevel > 0
    out.surfCorr = WW;
end

ball = decomposeBallLens(WW);% false triggers the correction of only the linear rotation part
MCorr = permute(ball.Mcorr,[1,3,2]);
MCorrT = bsxfun(@times,MCorr([1,4,7,2,5,8,3,6,9],:,:),[1;1;-1;1;1;-1;-1;-1;1]);
MM = MatrixMultiply(repmat(MCorrT,[1,size(MM,2),1]),MatrixMultiply(MM,repmat(MCorr,[1,size(MM,2),1])));

if outputLevel > 0
    out.ball = ball;
end
if outputLevel > 1
    out.OA3 = decomposeRot(MM);
end

%% now, extract signal of the outher sheath interface
% rrz = min(cath(3,:))-9:max(cath(3,:)) + 9;
% 
% ww = shiftdim(dop(rrz,:),-1);
% Mm = squeeze(bsxfun(@rdivide,sum(bsxfun(@times,MM(:,rrz,:),ww),2),sum(ww*9,2)));
% thcorr = decomposeRot(euclideanRotation(imfilter(Mm,ones(1,fw)/fw,'circular')));
% sheathRot = atan2(thcorr(2,:),thcorr(1,:));
% sheathRet = sqrt(sum(thcorr.^2,1));
% 
% if outputLevel > 0
%     out.sheath = sheathRot(:);
%     out.sheathRet = sheathRet(:);
% end
% 

% extract surface signal
mask = dop>0.8;
mask(1,:) = false; % make sure, at least the very first line is zero
%mask(1:refInds(3),:) = false; % make sure, at least the very first line is zero
mask(end,:) = false; % as well as the last one, to have an even pair of up and downward edges in each A-line

flipmask = flip(mask,1);
temp = cumsum(flipmask,1);
dm = cat(1,zeros(1,size(mask,2)),diff(flipmask,[],1));
lInds1 = find(dm>0);
lInds2 = find(dm<0);
ss = zeros(size(temp));
ss(lInds2) = temp(lInds2)-temp(lInds1);
[~,mp] = max(flip(ss,1),[],1);

% improve the surface detection by medfiltering and a 'smart' approach to
% pick the original and filtered surface signal
mpf = medfilt1(cat(2,mp,mp),fw);
stdmp = sqrt(medfilt1((cat(2,mp,mp)-mpf).^2,fw));
stdmp = circshift(stdmp(dim(2)/2 + (1:dim(2))),dim(2)/2);
mpf = circshift(mpf(dim(2)/2 + (1:dim(2))),dim(2)/2);
mp(abs(mp-mpf)>2*stdmp) = mpf(abs(mp-mpf)>2*stdmp);
mp = round(mp);
mp = max(mp + 8,cath(3,:));
%mp(mp<refInds(3)) = refInds(3);

inds = sub2ind(size(MM),ones(1,numel(mp)),mp,1:size(MM,3));
inds = bsxfun(@plus,inds,(0:size(MM,1)-1)');
Mm = MM(inds);
OAtissue = decomposeRot(euclideanRotation(imfilter(Mm,ones(1,fw)/fw,'circular')));

out.tissueAngle = atan2(OAtissue(2,:),OAtissue(1,:));
out.tissueRet = sqrt(sum(OAtissue.^2,1));
out.tissueSurf = mp;

[ang,err] = estimateOrientation(OAtissue,1);% mean 
out.sheathAngle = ang;
out.errSheathAngle = err;

% % use the signal of the sample surface to determine the best V-rotation
% % correction strategy
% [Vcorr,OAtissuep,ang,err,psi,type,ballErr,OAmodel] = refineVcorrection(ball,OAtissue);

%out.sheathAngle = ang;
%out.errSheathAngle = err;
%out.refineType = type;
%out.tissueOffset = psi;
%out.Vcorr = Vcorr;
%out.tissueAngle = squeeze(atan2(OAtissuep(2,:,:),OAtissuep(1,:,:)));
%out.tissueRet = squeeze(sqrt(sum(OAtissuep.^2,1)));
%out.ballErr = ballErr;
%out.OAmodel = OAmodel;

%Vcorr = Vcorr(2,:);% let's use the trace lsq solution
%%
% generate mask to exclude points with low DOP from correcting the local
% retardation
mmask = (dop>dopTh)&(bsxfun(@minus,(1:size(dop,1))',cath(3,:))>0);
%mmask(1:refInds(1)-1,:) = 0;

mmask = imdilate(imerode(mmask,ones(fwz,1)),ones(fwz,1));

dim = size(MM);
w = zeros([3,dim(2),dim(3)]);
if outputLevel>0
    w2 = zeros([3,dim(2),dim(3)]);
end

N = repmat([1;0;0;0;1;0;0;0;1],[1,1,dim(3)]);
%N = permute(makeRot(cat(1,zeros(2,numel(Vcorr)),-Vcorr)),[1,3,2]);

W = decomposeRot(MM);
Msqinv = makeRot(-W/2);
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
%        nloc = MatrixMultiply(MM(:,indz,:,:),MM([1,4,7,2,5,8,3,6,9],indz-1,:,:));%
%        %this results in a optic axis with a V-component!

        Wloc = decomposeRot(nloc);
        w2(:,indz,:) = squeeze(Wloc);
    end
end

if outputLevel>0
    % without depth correction
    Omegaf = imfilter(permute(w2,[2,3,1]),ones(fwz,1)/fwz);
    phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
    out.phi2 = phi2;
end
    
%Omega = bsxfun(@times,permute(pa,[2,3,1]),ret);
if ~cumulative % depth-resolved signal
    Omegaf = imfilter(permute(w,[2,3,1]),ones(fwz,1)/fwz);
    retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
    phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));

    Omegafn = bsxfun(@rdivide,Omegaf,sqrt(sum(Omegaf.^2,3)));
else % outputting the cumulative signal
    Omegaf = permute(W,[2,3,1]);
    retFinal = sqrt(sum(Omegaf.^2,3))/pi*100;% scale pi to 100
    phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));

    Omegafn = bsxfun(@rdivide,Omegaf,sqrt(sum(Omegaf.^2,3)));
end

out.dop = dop;
out.mask = mmask;
out.ret = retFinal;
out.phi = mod(phi + pi,2*pi)-pi;%do not include sheathAngle, as it may vary in between B-scans
%out.phi = mod(phi + out.sheathAngle + pi,2*pi)-pi;
out.PA = Omegafn;


out.fwz = fwz;
out.fwx = fwx;
out.int = int;
out.fw = fw;

