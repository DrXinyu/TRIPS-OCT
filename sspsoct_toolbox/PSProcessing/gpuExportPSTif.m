function TT = gpuExportPSTif(path,fout,N,fwx,dz,logLim,maxRet,opt)
% GPU version of exporting the pc reconstructed data to a .pstif file;
% uses the matlab internal GPU support
%
% For optimum performance, the B-scans are reconstructed for "splitSize"
% A-lines at a time.
%
% Arguments:
% path: points to the measurement folder
% fout: name of the output .pstif stack
% N: number of spectral bins to apply
% fwx: lateral filter fwhm
% dz: half width of axial offset for derivation of the local retardation
% logLim: lower and upper log limit for scaling of output tif
% maxRet: maximum local retardation value in deg/100um for scaling of
% output
% opt: optional structure with field .sliceInd to give a vector of the index
% slices to be reconstructed; otherwise all the available frames are
% reconstructed.
% opt: field .oop triggers out of plane averaging, by compounding the Stokes
% vectors from -oop:oop around the center slide.

% input parameters
rS.window = N;
ps.fwx = fwx;
ps.dz = dz;

% default settings
logF = readLogFile(path);
sliceInd = 1:logF.numImages;
oop = 0; % no oop averaging
If = false;

% optional settings
if nargin>7 && ~isempty(opt)
    if isstruct(opt)
        if isfield(opt,'sliceInd')
            sliceInd = opt.sliceInd;
        end
        if isfield(opt,'oop')
            oop = opt.oop;
        end
        if isfield(opt,'If')% output also the filtered intensity
            If = opt.If;
        end
    end
end

% clear output file
fid = fopen(fout,'w+');
fclose(fid);

logOS = logLim(1);
logD = diff(logLim);

% padd data already prior to call of gPSProcess for splitting of large
% frames
prepad = floor(ceil(1.5*fwx)/2)+floor(fwx/2);
postpad = ceil(ceil(1.5*fwx)/2)+ceil(fwx/2);
rS.prepad = prepad;
rS.postpad = postpad;
rS.logF = logF;

[map,disp] = readConfig(path);
bgr = readBgr(path);

rS.map = map;
rS.disp = disp;
rS.bgr = bgr;
st.map = map;
st.disp = disp;
st.bgr = bgr;

S1 = zeros([1024,logF.numAlinesImage/2,1,3]);
S2 = S1;


dim = size(S1);

splitSize = 512;
step = min(dim(2),splitSize);
if dim(2)>splitSize
    Split = ceil(dim(2)/splitSize);
else
    Split = 1;
end

if oop == 0
    %%% no oop averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tag to be written in the .pstif header to keep track of the
    % reconstruction parameters
    if ~If
        tag = sprintf('Format:\tIRD\nInt:\t%d\t%d\nRet:\t0\t%d\nDOP:\t0\t1\nParamsNFwxDz:\t%d\t%d\t%d\n%s',logLim(1),logLim(2),maxRet,N,fwx,dz,datestr(clock));
    else
        tag = sprintf('Format:\tIRDF\nInt:\t%d\t%d\nRet:\t0\t%d\nDOP:\t0\t1\nIntFilt:\t%d\t%d\nParamsNFwxDz:\t%d\t%d\t%d\n%s',logLim(1),logLim(2),maxRet,logLim(1),logLim(2),N,fwx,dz,datestr(clock));
    end
    
    for ind = 1:numel(sliceInd)

        % normal tomogram for intensity view
        tt1 = tic();
        st.gpu = true;
        st.logF = logF;
        [S1,S2] = recstrTom(path,[sliceInd(ind),sliceInd(ind)],st);
        ii = gather(tom2Int(S1,S2));

        ret = zeros([dim(1),step,Split]);
        dop = ret;
        iif = ret;
        
        % loop over splits
        for splitInd = 1:Split

            % pc Stokes vectors
            rS.gpu = true;
            rS.firstAline = 1 + (splitInd-1)*splitSize*2;
            rS.lastAline = splitSize*2+(splitInd-1)*splitSize*2;
            [gSpc1,gSpc2] = recstrTom(path,sliceInd(ind),rS);
            % ps reconstruction

            if ~If
                [rr,dd] = gpuPSProcess(gSpc1,gSpc2,fwx,dz,2*N-1,prepad,postpad);
            else
                [rr,dd,iiff] = gpuPSProcess(gSpc1,gSpc2,fwx,dz,2*N-1,prepad,postpad);
                iif(:,:,splitInd) = iiff;
            end
            
            ret(:,:,splitInd) = rr;
            dop(:,:,splitInd) = dd;
        end


        iiout = uint8(255*(10*log10(ii)-logOS)/logD);
        rrout = uint8(255*ret(:,:)/maxRet);
        ddout = uint8(255*dop(:,:));

        imwrite(iiout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
        imwrite(rrout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
        imwrite(ddout,fout,'tif','WriteMode','append','Compression','none','Description',tag);

        if If
            iifout = uint8(255*(10*log10(iif(:,:))-logOS)/logD);
            imwrite(iifout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
        end
        
        TT(ind) = toc(tt1)
        ind
    end

else
%%% with oop averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~If
        tag = sprintf('Format:\tIRD\nInt:\t%d\t%d\nRet:\t0\t%d\nDOP:\t0\t1\nParamsNFwxOopDz:\t%d\t%d\t%d\t%d\n%s',logLim(1),logLim(2),maxRet,N,fwx,2*oop+1,dz,datestr(clock));
    else
        tag = sprintf('Format:\tIRDF\nInt:\t%d\t%d\nRet:\t0\t%d\nDOP:\t0\t1\nIntFilt:\t%d\t%d\nParamsNFwxOopDz:\t%d\t%d\t%d\t%d\n%s',logLim(1),logLim(2),maxRet,logLim(1),logLim(2),N,fwx,2*oop+1,dz,datestr(clock));
    end
    % build up block of 2*oop + 1 Bscans
    blockInd = 0;
    for locSliceInd = max(sliceInd(1) + (-oop:oop),1)
        blockInd = blockInd + 1;
        for splitInd = 1:Split

            % pc Stokes vectors
            rS.gpu = true;
            rS.firstAline = 1 + (splitInd-1)*splitSize*2;
            rS.lastAline = splitSize*2+(splitInd-1)*splitSize*2;
            rS.fullStokes = true;
            [gSpc1loc,gSpc2loc] = recstrTom(path,[locSliceInd,locSliceInd],rS);
            gSpc1(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,blockInd) = gather(gSpc1loc);
            gSpc2(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,blockInd) = gather(gSpc2loc);
        end            
    end

    for ind = sliceInd
        tStart = tic();
        [S1,S2] = recstrTom(path,[ind,ind],st);
        ii = gather(tom2Int(S1,S2));
        dim = size(S1);


        ret = zeros([dim(1),step,Split]);
        dop = ret;
        iif = ret;

        % loop over splits
        for splitInd = 1:Split

            % pc Stokes vectors
            rS.gpu = true;
            rS.firstAline = 1 + (splitInd-1)*splitSize*2;
            rS.lastAline = splitSize*2+(splitInd-1)*splitSize*2;
            rS.fullStokes = true;
            [gSpc1loc,gSpc2loc] = recstrTom(path,min([ind,ind]+oop,logF.numImages),rS);
            gSpc1(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,blockInd) = gather(gSpc1loc);
            gSpc2(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,blockInd) = gather(gSpc2loc);
        end
        
        for splitInd = 1:Split

            if ~If
                [rr,dd] = gpuPSProcessOOP(gpuArray(mean(gSpc1(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,:),5)),gpuArray(mean(gSpc2(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,:),5)),fwx,dz,2*N-1,prepad,postpad);
            else
                [rr,dd,iiff] = gpuPSProcessOOP(gpuArray(mean(gSpc1(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,:),5)),gpuArray(mean(gSpc2(:,( 1 + (splitInd-1)*(splitSize)):(splitSize+rS.prepad+rS.postpad+(splitInd-1)*splitSize),:,:,:),5)),fwx,dz,2*N-1,prepad,postpad);
                iif(:,:,splitInd) = iiff;
            end
            
            ret(:,:,splitInd) = rr;
            dop(:,:,splitInd) = dd;
        end


        logOS = logLim(1);
        logD = diff(logLim);

        iiout = uint8(255*(10*log10(ii)-logOS)/logD);
        rrout = uint8(255*ret(:,:)/maxRet);
        ddout = uint8(255*dop(:,:));

        imwrite(iiout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
        imwrite(rrout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
        imwrite(ddout,fout,'tif','WriteMode','append','Compression','none','Description',tag);

        if If
            iifout = uint8(255*(10*log10(iif(:,:))-logOS)/logD);
            imwrite(iifout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
        end
        slideTime = toc(tStart)
        blockInd = mod(blockInd,2*oop+1)+1;
    end
end


    
    