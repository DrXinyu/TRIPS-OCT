function [Msqinv,Wout,mmnew] = unwrapOAx(MM,dop,dopTh,mmIn)
%Wout = unwrapOAx(W1,dop,dopTh) generates the most plausible unwrapped
%optic axis Wout, based on the original wrapped W1 and the dop signal,
%which serves to ignore areas with dop<dopTh.

if nargin<3 || isempty(dopTh)
    dopTh = .6;
end
if nargin<4 || isempty(mmIn)
    mmIn = [];
end
dopm = (dop>dopTh);

if size(MM,1)==3
    W1 = MM;
    MM = makeRot(W1);
else
    W1 = decomposeRot(MM);
end

Wn = W1./sqrt(sum(W1.^2,1));
W2 = W1 - 2*pi*Wn; 

%Compute differential matrix
dr = decomposeRot(MatrixMultiply(MM(:,2:end,:),MM([1;4;7;2;5;8;3;6;9],1:end-1,:)));
% comput scalar product between dr and dw
dw = (W1(:,2:end,:)-W1(:,1:end-1,:));
ss1 = sign(sum(dw(1:2,:,:).*dr(1:2,:,:),1) + sqrt(eps))<0;% add tiny number to avoid rounding issues;third dimension has to be zero
ss1(sum(dw.^2,1)>6.1^2) = 1; % force points with dw close to 2*pi to wrap.

% initial unwrapping mask along axial direction
mm = squeeze(mod(cumsum(cat(2,zeros(1,1,size(W1,3)),ss1),2),2)>0);

% unwrapping mask in lateral direction
dr = decomposeRot(MatrixMultiply(MM(:,:,2:end),MM([1;4;7;2;5;8;3;6;9],:,1:end-1)));
% comput scalar product between dr and dw
dw = (W1(:,:,2:end)-W1(:,:,1:end-1));
ss2 = sign(sum(dw(1:2,:,:).*dr(1:2,:,:),1) + sqrt(eps))<0;% add tiny number to avoid rounding issues
ss2(sum(dw.^2,1)>6.1^2) = 1; % force points with dw close to 2*pi to wrap.
drn2 = squeeze(sqrt(sum(dr.^2,1)));


edgeMap = cat(1,dopm(1,:),diff(dopm,[],1)==1);
edgeEndMap = cat(1,diff(dopm,[],1)==-1,dopm(end,:));

% index each segment and determine its length
temp = cumsum(dopm);
inds = find(edgeMap);
ll = temp(edgeEndMap) - temp(edgeMap) + 1;

% generate label map (each line segment of dop>dopTh gets a label)
lb = NaN(size(dop));
for ind = 1:numel(inds)
    lb(inds(ind) + (0:ll(ind)-1)) = ind;
end

mmxor = xor(xor(mm(:,1:end-1),mm(:,2:end)),squeeze(ss2));% wherever mmxor there is an inconsistency in the current mm
errApp = drn2.*(~mmxor) + (2*pi-drn2).*mmxor - pi; %this is in the range of
%-pi to pi and designates the 'distance' between the two rotation matrices,
%with a negative value indicating correct wrapping and a positive value a
%inconsistent wrapping. (0 to 2*pi is more logical, but centering it on
%zeros allows taking the negative value of sums to get the wrapped-value).

signs = NaN(1,numel(inds));% indicates the requirement to unwrap the current A-line segment (1) or not (0).
seqLabel = zeros(1,numel(inds));% assigns label of contiguous internally consistently unwrapped areas to each A-line segment
Nz = size(dop,1);
Nlines = size(dop,2);
totSegInd = 0;
conflictInds = [];
for ind = 1:numel(inds)
    if isnan(signs(ind))
        totSegInd = totSegInd + 1;
        signs(ind) = 0;
        seqLabel(ind) = totSegInd; 
        try
            [seqOut,status,maxError,maxErrorInd] = sequentialUnwrapping(ind,[],ind,1e6,NaN);
        catch
            disp('asdf');
            seqOut
        end
        labelSize(totSegInd) = sum(ll(seqLabel == totSegInd));
    end
end

% generate new unwrapping mask
mmnew = zeros(size(mm));
labelMap = zeros(size(mm));
for ind = 1:numel(inds)
    mmnew(inds(ind) + (0:ll(ind)-1)) = mod((mm(inds(ind) + (0:ll(ind)-1)))+signs(ind),2);
    labelMap(inds(ind) + (0:ll(ind)-1)) = seqLabel(ind);
end

mmxor = xor(xor(mmnew(:,1:end-1),mmnew(:,2:end)),squeeze(ss2));% wherever mmxor there is an error in the current mm
errAppFinal = drn2.*(~mmxor) + (2*pi-drn2).*mmxor; %Here, the range 0 to 2*pi was used


% merge the individual sequences
% simply repeat closest W1 in areas of dop<dopTh
lastLinelM = NaN(1,size(dopm,2));
lastLinemm = NaN(1,size(dopm,2));
lMp = labelMap;
mmp = mmnew;
for ind = 1:size(dopm,1)
    lMp(ind,~dopm(ind,:)) = lastLinelM(1,~dopm(ind,:));
    lastLinelM = lMp(ind,:);
    mmp(ind,~dopm(ind,:)) = lastLinemm(1,~dopm(ind,:));
    lastLinemm = mmp(ind,:);
end

% find segment interfaces
interFaces = lMp.*cat(1,abs(diff(lMp,[],1))>0,false(1,size(dopm,2)));

% check out border of each label
HTable = zeros(totSegInd);
for ind = 1:totSegInd
    inds = find(interFaces==ind);
    % identify neighboring segments
    nbs = unique(lMp(inds+1));
    for jnd = 1:numel(nbs)
        jnds = inds(lMp(inds+1) == nbs(jnd));
        HTable(nbs(jnd),ind) =  sum(.5-xor(mmp(jnds),mmp(jnds+1)));
    end
end
% remove loops from table, where two segments depend on each other
% (enclosed islands);
[inds,jnds] = ind2sub(size(HTable),find(((abs(HTable)>0) + (abs(HTable)>0)')==2));
if ~isempty(inds)
    HTable(sub2ind(size(HTable),inds(inds<jnds),jnds(inds<jnds))) = 0;
end

% get sign-dependency for all segments
ss = findSignDependency(HTable);
sstot = sum(ss,1);% if ss has multiple rows, there are non-connected areas in the tomogram; for now, let's ignore this case.


ss = sum(sstot,1);
for ind = 1:numel(ss)
    if ss(ind) == -1
        mmnew = mod(mmnew + (labelMap == ind),2);
    end
end

if ~isempty(mmIn)
    if sum(sum(xor(~mmnew,mmIn).*dopm)) < sum(sum(xor(mmnew,mmIn).*dopm))
        mmnew = ~mmnew;
    end
end

Wout = W1;
Wout(1,mmnew>0) = W2(1,mmnew>0);
Wout(2,mmnew>0) = W2(2,mmnew>0);
Wout(3,mmnew>0) = W2(3,mmnew>0);

% What we really need is the square root of the input rotation matrix
Msqinv = makeRot(-Wout/2);


function [seqOut,status,maxError,maxErrorInd] = sequentialUnwrapping(ind,callerInd,seqIn,maxErrorIn,maxErrorIndIn)
%[seqOut,status,maxError,maxErrorInd] =
%sequentialUnwrapping(ind,callerInd,seqIn,maxErrorIn,maxErrorIndIn) finds
%all neighboring A-line segments of the segment labelled 'ind' in the label
%map 'lb'. 'callerInd' is the calling segment, which is also a neighbor.
%'seqIn' is the already established sequence of segment indices in the
%current 'branch'. 'maxErrorIn' and 'maxErrorIndIn' indicate the element in
%the current sequence that resulted in the largest error, and is the
%preferred 'breakingPoint' within that segment if there is a conflict
%arising (i.e. looped back sequence does not match in sign).

%seqOut: list of points in the current sequence
%status: Indicates status of current (growing) edge:
% status = 0 : reached end point (no new neighbors)
% status = 1 : reached existing node
% status = -1 : reached existing node with conflict
% status = 2 : reached new node
%
%maxError: largest error along the current sequence
%maxErrorInd: corresponding segment index

% identify all neighbors
if inds(ind)<(Nlines-1)*Nz
    rightNeighbors = min(lb((inds(ind) + (0:ll(ind)-1))+Nz)):max(lb((inds(ind) + (0:ll(ind)-1))+Nz));
else
    rightNeighbors = [];
end
if inds(ind)>Nz
    leftNeighbors = min(lb((inds(ind) + (0:ll(ind)-1))-Nz)):max(lb((inds(ind) + (0:ll(ind)-1))-Nz));
else
    leftNeighbors = [];
end

% remove NaNs
rightNeighbors = rightNeighbors(~isnan(rightNeighbors));
leftNeighbors = leftNeighbors(~isnan(leftNeighbors));
neighbors = cat(1,rightNeighbors(:),leftNeighbors(:));
offsetBool = cat(1,true(numel(rightNeighbors),1),false(numel(leftNeighbors),1));% needed to efficiently assess left and right hand side neighbors

% remove the callerInd to only have true neighbors
[neighbors,ids] = setdiff(neighbors,callerInd);
offsetBool = offsetBool(ids);

if numel(neighbors)==0 && numel(seqIn)==1% special case; segment that has no neighbors; return single node
%    segSeq = cat(2,segSeq,ind);
    seqOut = seqIn;
    status = NaN;
    maxError = maxErrorIn;
    maxErrorInd = maxErrorIndIn;
elseif numel(neighbors)==0% no more neighbors; endPoint
    seqOut = seqIn;
    status = 0;% reached end point
    maxError = maxErrorIn; % no update
    maxErrorInd = maxErrorIndIn; % no update
elseif numel(neighbors) > 1 || numel(seqIn)==1% new node, iterate over all neighbors
    for jnd = 1:numel(neighbors)
        if isnan(signs(neighbors(jnd)))% only go down this branch if it has not already been 'signed'
            startInd = max(inds(ind)-Nz*~offsetBool(jnd),inds(neighbors(jnd))-Nz*offsetBool(jnd));
            endInd = max(inds(ind) + ll(ind)-1 - Nz*~offsetBool(jnd),inds(neighbors(jnd))-Nz*offsetBool(jnd) + ll(neighbors(jnd))-1);
            errSum = sum(errApp(startInd:endInd));

            maxError = errSum;
            maxErrorInd = neighbors(jnd);
            
            if errSum<=0
                signs(neighbors(jnd)) = signs(ind);
            else
                signs(neighbors(jnd)) = mod(signs(ind)+1,2);
                maxError = -errSum;
            end
            seqLabel(neighbors(jnd)) = seqLabel(ind);

            % recursive call
            [seqOut,status,maxErrorBranch,maxErrorIndBranch] = sequentialUnwrapping(neighbors(jnd),ind,[ind,neighbors(jnd)],maxError,maxErrorInd);
            % make node table; for each neighbor, report the connecting
            if status == -1 
                startInd = find(seqOut == maxErrorIndBranch);
                for locInd = startInd:numel(seqOut)
                    signs(seqOut(locInd)) = mod(signs(seqOut(locInd))+1,2);
                end
                conflictInds = cat(2,conflictInds,maxErrorIndBranch);    
            end
            % node, the lowest cost and its index
%            segSeq = cat(2,segSeq,seqOut);

            % report back 
            seqOut = seqIn;%cat(2,seqIn,neighbors(jnd));
            status = 2;% reached new node
            maxError = maxErrorIn; % no update
            maxErrorInd = maxErrorIndIn; % no update
        end
    end
elseif numel(neighbors)==1% one new neighbor; stop if looping back, otherwise iterate
    jnd = 1;
    startInd = max(inds(ind)-Nz*~offsetBool(jnd),inds(neighbors(jnd))-Nz*offsetBool(jnd));
    endInd = max(inds(ind) + ll(ind)-1 - Nz*~offsetBool(jnd),inds(neighbors(jnd))-Nz*offsetBool(jnd) + ll(neighbors(jnd))-1);
    errSum = sum(errApp(startInd:endInd));

    if min(errSum,-errSum)>maxErrorIn
        maxError = min(errSum,-errSum);
        maxErrorInd = neighbors(jnd);
    else
        maxError = maxErrorIn;
        maxErrorInd = maxErrorIndIn;
    end

    seqOut = cat(2,seqIn,neighbors(jnd));
    if isnan(signs(neighbors(jnd)))% not yet 'signed'; iterate
        if errSum<=0
            signs(neighbors(jnd)) = signs(ind);
        else
            signs(neighbors(jnd)) = mod(signs(ind)+1,2);
        end
        seqLabel(neighbors(jnd)) = seqLabel(ind);
        % recursive call
        [seqOut,status,maxError,maxErrorInd] = sequentialUnwrapping(neighbors(jnd),ind,seqOut,maxError,maxErrorInd);
    else % looping back
        seqOut = seqIn;% segment is already in sequence
        if xor(errSum<=0,xor(signs(ind), signs(neighbors(jnd))))
            status = 1;%neighbors(jnd);
        else
            status = -1;%-neighbors(jnd);
        end
    end
end

end% function sequentialUnwrapping

end% function unwrapOAX


