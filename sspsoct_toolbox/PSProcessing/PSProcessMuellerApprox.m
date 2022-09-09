function out= PSProcessMuellerApprox(tom,procStruct)
%out= PSProcessMuellerApprox(tom,procStruct) reconstructs the local
%retardation and depolarization of the freqmux tomogram in tom, where tom
%has the four Jones compoments [tom.ch11,tom.ch12;tom.ch21,tomch22]. THe
%processing is performed following M. Villiger, et al., "Deep tissue volume
%imaging of birefringence through fibre-optic needle probes for the 
%delineation of breast tumour," Sci. Rep. 6, 28771 (2016).
%Also runs when called with gpuArray initalized variables.
%
%Modification on Jan 3, 2019: Also return the conventional DOP, when
%procStruct.DOP is set to true

Mull = Jones2MuellerTom(cat(3,tom.ch11,tom.ch12,tom.ch21,tom.ch22));

% default parameters
dz = 4;
fwx = 6;
fwz = 4;
dzres = 4.3;%um/px

% input parsing
if isfield(procStruct,'fwz')
    fwz = procStruct.fwz;
end

if isfield(procStruct,'fwx')
    fwx = procStruct.fwx;
end

if isfield(procStruct,'dz')
    dz = procStruct.dz;
end

if isfield(procStruct,'dzres')
    dzres = procStruct.dzres;
end


Nz = size(Mull,1);
NAlines = size(Mull,2);

% construct filter
nx = (round(fwx*1.5)-1)/2;
nz = (round(fwz*1.5)-1)/2;
nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
nz = linspace(-nz,nz,round(fwz*1.5))'*2*sqrt(log(2))/fwz;
h = exp(-nz.^2)*exp(-nx.^2);
h = h/sum(h(:));

Mullf = imfilter(Mull,h,'replicate')*1e-8;

% if requested, compute conventional DOP from S1 and S2, which correspond
% to Mullf*[1;1;0;0] and Mullf*[1;-1;0;0], respectively
if isfield(procStruct,'DOP') && procStruct.DOP
    out.DOP1 = sqrt(sum((Mullf(:,:,2:4) + Mullf(:,:,6:8)).^2,3))./(Mullf(:,:,1) + Mullf(:,:,5));
    out.DOP2 = sqrt(sum((Mullf(:,:,2:4) - Mullf(:,:,6:8)).^2,3))./(Mullf(:,:,1) - Mullf(:,:,5));
end
if isfield(procStruct,'If') && procStruct.If
    out.If = Mullf(:,:,1)*1e8;
end

if isfield(procStruct,'decomp') && procStruct.decomp
    %first, find diattenuation
    decomp = fastPolarDecomp(permute(Mullf,[3,1,2]));
    out.decomp = decomp;
end

% normalize Mullf by its determinant
posInds = [6,11,16,7,12,14,8,10,15; 5,11,16,7,12,13,8,9,15; 5,10,16,6,12,13,8,9,14; 5,10,15,6,11,13,7,9,14];
negInds = [8,11,14,6,12,15,7,10,16; 8,11,13,5,12,15,7,9,16; 8,10,13,5,12,14,6,9,16; 7,10,13,5,11,14,6,9,15];
preInds = [1,2,3,4]; 
ssign = [1,-1,1,-1];

DD = zeros(Nz,NAlines);
for ind = 1:4
    p = Mullf(:,:,posInds(ind,1)).*Mullf(:,:,posInds(ind,2)).*Mullf(:,:,posInds(ind,3))...
    + Mullf(:,:,posInds(ind,4)).*Mullf(:,:,posInds(ind,5)).*Mullf(:,:,posInds(ind,6))...
    + Mullf(:,:,posInds(ind,7)).*Mullf(:,:,posInds(ind,8)).*Mullf(:,:,posInds(ind,9));

    n = Mullf(:,:,negInds(ind,1)).*Mullf(:,:,negInds(ind,2)).*Mullf(:,:,negInds(ind,3))...
    + Mullf(:,:,negInds(ind,4)).*Mullf(:,:,negInds(ind,5)).*Mullf(:,:,negInds(ind,6))...
    + Mullf(:,:,negInds(ind,7)).*Mullf(:,:,negInds(ind,8)).*Mullf(:,:,negInds(ind,9));

    DD = DD + ssign(ind)*Mullf(:,:,preInds(ind)).*(p-n);
end

Mullf = bsxfun(@rdivide,Mullf,(abs(DD).^(1/4)));
Dep = sqrt((sum(Mullf.^2,3) - Mullf(:,:,1).^2)./(Mullf(:,:,1).^2)/3);

M1 = Mullf;%circshift(Mullf,Param.dz);
M2 = circshift(Mullf,-1);

% construct inverse of M1
M1inv = M1(:,:,[1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16]);
M1inv(:,:,[2,3,4,5,9,13]) = -M1inv(:,:,[2,3,4,5,9,13]);

for ind = 1:4
    for jnd = 1:4
        MM(:,:,ind + (jnd-1)*4) = sum(M2(:,:,ind + [0,4,8,12]).*M1inv(:,:,(jnd-1)*4 + [1,2,3,4]),3);
    end
end

Ret(:,:,1) = (MM(:,:,12)-MM(:,:,15))/2; 
Ret(:,:,2) = (MM(:,:,14)-MM(:,:,8))/2; 
Ret(:,:,3) = (MM(:,:,7)-MM(:,:,10))/2; 


Diat(:,:,1) = (MM(:,:,2)-MM(:,:,5))/2; 
Diat(:,:,2) = (MM(:,:,3)-MM(:,:,9))/2; 
Diat(:,:,3) = (MM(:,:,4)-MM(:,:,13))/2; 


hh = ones(dz*2,1);
% correcting for influence of Dep
Ret = imfilter(Ret,hh);
Diat = imfilter(Diat,hh);
%Ret = imfilter(bsxfun(@rdivide,Ret,sqrt(Dep)),hh);
%Diat = imfilter(bsxfun(@rdivide,Diat,sqrt(Dep)),hh);

out.ret = sqrt(sum(Ret.^2,3))/2/dz/dzres/pi*180*100;
out.diatt = sqrt(sum(Diat.^2,3))/2/dz/dzres;
out.dep = Dep;
out.PA = bsxfun(@rdivide,Ret,sqrt(sum(Ret.^2,3)));

