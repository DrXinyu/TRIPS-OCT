function out = exportOAPSTifNew(path,fout,fwx,dz,logLim,maxRet,N,opt,boolCatheter)
% script to generate OA pstifs using the newer reconstruction strategy
% without simple iterative correction.
% Under development!

if nargin<4 || isempty(dz)
    dz = 5;
end
if nargin<5 || isempty(logLim)
    logLim = [55,110];
end
if nargin<6 || isempty(maxRet)
    maxRet = 100;% degrees per 100µm
end
if nargin<7 || isempty(N)
    N = 5;
end
if nargin<9 || isempty(boolCatheter)
    boolCatheter = false;
end

% default settings
logF = readLogFile(path);
sliceInd = 1:logF.numImages;
st = struct;
[m,d] = readConfig(path);
dzres = 4.33;
boolphi2 = false;
fwaxial = 1; % exceptional axial filtering in absence of spectral binning
cumulative = false;% flag for processing cumulative signal, if requested
fwy = .5; % filtering in oop direction; 0.5 corresponds to a single slice

% optional settings
if nargin>7 && ~isempty(opt)
    if isstruct(opt)
        if isfield(opt,'sliceInd')
            sliceInd = opt.sliceInd;
        end
        if isfield(opt,'disp')
            st.disp = d.*opt.disp;
        end
        if isfield(opt,'dzres')
            dzres = st.dzres;
        end
        if isfield(opt,'phi2')
            boolphi2 = true;
        end
        if isfield(opt,'fwaxial')
            fwaxial = opt.fwaxial;
        end
        if isfield(opt,'cumulative')
            cumulative = opt.cumulative;
        end
        if isfield(opt,'fwy')% out of plane averaging
            fwy = opt.fwy;
        end
    end
end

if cumulative
    boolphi2 = true;
end
% construction of filter for out of plane averaging
ny = (round(fwy*1.5)-1)/2;
ny = linspace(-ny,ny,round(fwy*1.5))*2*sqrt(log(2))/fwy;
hy = exp(-ny.^2);
hy = shiftdim(hy/sum(hy(:)),-1);
oop = numel(hy);% number of slices required for out of plane averaging


fid = fopen(fout,'w+');
fclose(fid);

% construction of the tag for the header of the tif file; special,
% including two axis images
if boolphi2
    format = 'IRDPP';
    if ~cumulative
        if oop == 1
            Format = {'Int:',logLim(1),logLim(2),[];'Ret:',0,maxRet,[];'DOP:',0,1,[];'Phi:',-pi,pi,[];'Phi:',-pi,pi,[];'ParamsNFwxDz',N,fwx,dz};
        else
            Format = {'Int:',logLim(1),logLim(2),[],[];'Ret:',0,maxRet,[],[];'DOP:',0,1,[],[];'Phi:',-pi,pi,[],[];'Phi:',-pi,pi,[],[];'ParamsNFwxFwyDz',N,fwx,fwy,dz};
        end
    else
        if oop == 1
            Format = {'Int:',logLim(1),logLim(2),[];'Ret:',0,maxRet,[];'DOP:',0,1,[];'Phi:',-pi,pi,[];'Cum:',-pi,pi,[];'ParamsNFwxDz',N,fwx,dz};
        else
            Format = {'Int:',logLim(1),logLim(2),[],[];'Ret:',0,maxRet,[],[];'DOP:',0,1,[],[];'Phi:',-pi,pi,[],[];'Cum:',-pi,pi,[],[];'ParamsNFwxFwyDz',N,fwx,fwy,dz};
        end
    end
else
    format = 'IRDP';
    if oop == 1
        Format = {'Int:',logLim(1),logLim(2),[];'Ret:',0,maxRet,[];'DOP:',0,1,[];'Phi:',-pi,pi,[];'ParamsNFwxDz',N,fwx,dz};
    else
        Format = {'Int:',logLim(1),logLim(2),[],[];'Ret:',0,maxRet,[],[];'DOP:',0,1,[],[];'Phi:',-pi,pi,[],[];'ParamsNFwxFwyDz',N,fwx,fwy,dz};
    end
end
%    tt.Format = {'Int:',logLims(1),logLims(2),[];'Ret:',0,maxRet,[];'DOP:',0,1,[];'Phi:',0,pi,[];'Phi:',0,pi,[];'ParamsNFwxDz',5,12,5};
tag = sprintf('Format:\t%s\n',format);
for ind = 1:numel(format)
    tag = [tag,sprintf('%s\t%d\t%d\n',Format{ind,1},Format{ind,2},Format{ind,3})];
end
if size(Format,1)>numel(format)
    tag = [tag,sprintf('%s',Format{end,1})];
    for ind = 2:size(Format,2)
        tag = [tag,sprintf('\t%d',Format{end,ind})];
    end
    tag = [tag,sprintf('\n%s',datestr(clock))];
end

out = struct;
SS1w = [];
SS2w = [];
lastInds = [];
for ind = 1:numel(sliceInd)
    
    t = tic;

    [S1,S2] = recstrTom(path,[1,1]*sliceInd(ind),st);
    int = tom2Int(S1,S2);

    % spectral binning reconstruction
    st2 = st;
    st2.window = N;
    st2.skipLastPoints = 30;
    st2.fullStokes = true;

    if oop==1
        [S1w,S2w] = recstrTom(path,[1,1]*sliceInd(ind),st2);
    else
        % these are the slice indices that we need
        inds = ind + (-floor(oop/2):floor((oop-1)/2));
        indsmask = inds>0&inds<=numel(sliceInd);% they have to be valid
        inds = setdiff(sliceInd(inds(indsmask)),lastInds);% difference to previously reconstructed slices
        if ~isempty(inds)
            [s1w,s2w] = recstrTom(path,inds,st2);
            % permute SS1w if inds was a scalar
            order = 1:numel(size(s1w));
            if numel(inds) == 1
                order = [order(1:2),numel(order)+1,order(3:end)];
            end
            % figure out slices to recycle
            startInd = max((size(SS1w,3)-oop + numel(inds)+1),1);
            SS1w = cat(3,SS1w(:,:,startInd:end,:,:),permute(s1w,order));
            SS2w = cat(3,SS2w(:,:,startInd:end,:,:),permute(s2w,order));
        else
            SS1w = SS1w(:,:,2:end,:,:);
            SS2w = SS2w(:,:,2:end,:,:);
        end       
        S1w = squeeze(sum(SS1w.*hy(:,:,indsmask),3));
        S2w = squeeze(sum(SS2w.*hy(:,:,indsmask),3));
        lastInds = cat(2,lastInds,inds);
    end
    
    pstruct = struct;
    pstruct.fwx = fwx;
    pstruct.fwz = dz;
    pstruct.dzres = dzres;
    pstruct.fwaxial = fwaxial;
    pstruct.cumulative = cumulative;
    if isfield(out,'wcorrEff') % if this is a further iteration, recycle the previously determined wcorr
        pstruct.wcorr = out.wcorrEff;
    end
    if isfield(out,'rcEff') % if this is a further iteration, recycle the previously determined wcorr
        pstruct.rc = out.rcEff;
    end
    if isfield(out,'unwrapMask') % try to ensure continuity of unwrapping
        pstruct.unwrapMask = out.unwrapMask;
    end
    
    out = OAProcessing(S1w,S2w,pstruct,1);

%     if boolCatheter
%     else
%     end
    
    byteString = writeToUint8(out);
%    oaStr = readFromUint8(byteString,9,512)
    
    intout = uint8(255*(10*log10(int)-logLim(1))/diff(logLim));
    retout = uint8(255*out.ret/maxRet);
    dopout = uint8(255*out.dop);
    phiout = uint8(255*(out.phi+pi)/2/pi)';
    % add byteString to phiout
    phiout(1:numel(byteString)) = byteString;
    phiout = phiout';

    imwrite(intout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    imwrite(retout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    imwrite(dopout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    imwrite(phiout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    if boolphi2
        phi2out = uint8(255*(out.phi2+pi)/2/pi);
        imwrite(phi2out,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    end

    dtime = toc(t);
    
    sprintf('Section %d processed in %.2f s\n',sliceInd(ind),dtime)
end
