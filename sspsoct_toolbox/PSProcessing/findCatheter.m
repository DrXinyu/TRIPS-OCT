function lines = findCatheter(tom,Npoints)
% lines = findCatheter(tom) finds the position of the ball lens and the
% catheter interfaces and returns them in lines.
% It uses an edge detection filter to identify candidate pixels and then
% combines edge points that roughly line up to identify lines that
% continuously go across the entire B-scan.
% There are some additional criteria that enable identification of the
% correct lines from the candidate lines.

if nargin<2
    Npoints = 50;
end

maxZPos = 512;
maxCathWidth = 35;

% construct a filter that detects edges
fwz = 25;
fwx = 20;
sig1 = .1;
sig2 = 1/3;
ww = .5;
ths = .2;

hh1 = repmat((1+ww)*exp(-linspace(-1,1,fwz)'.^2/sig1^2)-ww*exp(-linspace(-1,1,fwz)'.^2/sig2^2),1,fwx);
hh1 = hh1/sqrt(sum(sum(hh1.^2)));

hh2 = ones(fwz,fwx);
gg1 = imfilter(tom,hh1,'circular');
gg2 = sqrt(imfilter(tom.^2,hh2,'circular'));
temp = gg1./gg2;

mask1 = temp>ths;
mask2 = cat(1,diff(sign(cat(1,zeros(1,size(temp,2)),diff(temp,[],1))),[],1),zeros(1,size(temp,2)))==-2;
mask = mask1.*mask2;

xlocs = round(linspace(1,size(tom,2),Npoints));
inds = find(mask(:,xlocs)==1);
[PosZ,xinds] = ind2sub([size(tom,1),numel(xlocs)],inds);
PosX = xlocs(xinds);
PosX = PosX(PosZ<maxZPos);
PosZ = PosZ(PosZ<maxZPos);

% find right hand neighbours
maxAxOffset = 4;
for pInd = 1:numel(PosZ)
    pp = find(abs(PosZ(pInd+1:end)-PosZ(pInd))<=maxAxOffset,1,'first');
    if ~isempty(pp)
        rh(pInd) = pp + pInd;
    else
        rh(pInd) = NaN;
    end
end

% left neighbours
for pInd = 1:numel(PosZ)
    pp = find(abs(PosZ(1:pInd-1)-PosZ(pInd))<=maxAxOffset,1,'last');
    if ~isempty(pp)
        lh(pInd) = pp;
    else
        lh(pInd) = NaN;
    end
end

% find way points, including those with NaN on the left
wp = [];
for pInd = 1:numel(PosZ)
    if isnan(lh(pInd)) || rh(lh(pInd))~=pInd
        wp = cat(1,wp,pInd);
    end
end

pEl = [];
for pInd = 1:numel(wp)

    % first go to right
    pnext = wp(pInd);
    Plist = pnext;
    while ~isnan(rh(pnext))
        pnext = rh(pnext);
        Plist = cat(2,Plist,pnext);
    end
    
    pnext = wp(pInd);
    while ~isnan(lh(pnext))
        pnext = lh(pnext);
        Plist = cat(2,pnext,Plist);
    end
    
    PList{pInd} = Plist;
    NN(pInd) = numel(Plist);
    WW(pInd) = PosX(Plist(end))-PosX(Plist(1));
    ZZ(pInd) = mean(PosZ(Plist(end)));
    
    StartPoint(pInd) = Plist(1);
    EndPoint(pInd) = Plist(end);
    
    %check if this starting and end points have been used previously;
    %eliminate poorer line
    pp = find(StartPoint(1:end-1)==StartPoint(end));
    if ~isempty(pp)
        
        if WW(pInd)>WW(pp) || std(diff(PosZ(Plist)))<std(diff(PosZ(PList{pp})))
            % pInd is better, eliminate pp
            pEl = cat(2,pEl,pp);
        else
            % pp is better, eliminate pInd
            pEl = cat(2,pEl,pInd);
        end
        StartPoint(pEl(end)) = NaN;
        EndPoint(pEl(end)) = NaN;
    end
    pp = find(EndPoint(1:end-1)==EndPoint(end));
    if ~isempty(pp)
        
        if WW(pInd)>WW(pp) || std(diff(PosZ(Plist)))<std(diff(PosZ(PList{pp})))
            % pInd is better, eliminate pp
            pEl = cat(2,pEl,pp);
        else
            % pp is better, eliminate pInd
            pEl = cat(2,pEl,pInd);
        end
        StartPoint(pEl(end)) = NaN;
        EndPoint(pEl(end)) = NaN;
    end
end

pinds = setdiff(1:numel(wp),pEl);
PList = PList(pinds);
NN = NN(pinds);
WW = WW(pinds);
ZZ = ZZ(pinds);

[~,pp] = sort(NN,'descend');
candInds = pp(1:max(find(NN(pp)>Npoints*.6,1,'last'),3));%  at least three points or more than 60% of Npoints
%candInds = find(NN>Npoints*.6);
[~,pp] = sort(ZZ(candInds),'ascend');
candInds = candInds(pp);
%candInds = candInds(pp(max(numel(candInds)-2,1):end));


for ind = 1:numel(candInds)
    lines(ind,:) = round(interp1(PosX(PList{candInds(ind)}),PosZ(PList{candInds(ind)}),1:size(tom,2),'linear','extrap'));
    
    lInds(ind,:) = sub2ind(size(tom),lines(ind,:),1:size(tom,2));

    peakSignal(ind) = median(tom(lInds(ind,:)));
end

for ind = 1:numel(candInds)-1
    midSignal(ind) = median(tom(lInds(ind,:) + round(min(mean(diff(lines(ind:ind+1,:),[],1))/2,maxCathWidth/2))));
%   midSignal(ind-1) = median(tom(lInds - round(mean(diff(lines(ind-1:ind,:),[],1))/2)));
end

if numel(candInds)>=3
    % find the true ball lens signal
    pp = find((midSignal(2:end)>midSignal(1:end-1))&(midSignal(2:end)>((peakSignal(2:end-1)+peakSignal(3:end))/2)/10^(25.2/10))&midSignal(2:end)>10^(70/10));
    
    if isempty(pp) % last line must be inner catheter sheath, so put pp to last but two
        pp = numel(candInds)-1;
    end
    
    lines = lines(pp:min(pp+2,numel(candInds)),:);
    candInds = candInds(pp:min(pp+2,numel(candInds)));
%     if midSignal(2)<midSignal(1) || midSignal(2)<mean(peakSignal(2:3))/10^(25/10) || midSignal(2)<10^(70/10)
%         candInds = candInds(2:3);
%         lines = lines(2:3,:);
%     end
end
% 
% % if only two lines detected, look for neigbouring points
% if numel(candInds) == 2
%     dz = PosZ(min(PList{candInds(2)}+1,numel(PosZ)))-PosZ(PList{candInds(2)});
%     dz = dz(dz>0);
%     lines(3,:) = lines(2,:) + median(dz);
% end

matchInd = PList{candInds(2)};
dz = PosZ(min(matchInd+1,numel(PosZ)))-PosZ(matchInd);
matchInd = matchInd((dz<maxCathWidth).*(dz>0)>0);
dz = dz((dz<maxCathWidth).*(dz>0)>0);
mdz = median(dz);
sdz = sqrt(median((dz-mdz).^2));
matchInd = matchInd(abs(dz-mdz)<4*max(sdz,1));
lines(3,:) = round(interp1(PosX(matchInd+1),PosZ(matchInd+1),1:size(tom,2),'linear','extrap'));
