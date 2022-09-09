function sstot = findSignDependency(HTable)
% HTable is a 'dependency' table showing the sign-relation between various
% indices. We assume that HTable has no circular dependence (i.e. element
% ij is zero if ji has a value. This is enforce by giving priority to lower
% triangular values.

HTableRed = HTable;
signTable = zeros(size(HTable));

% First, identify rows n with a single dependence on m. Add the n-th column
% to the m-th column (moving dependencies on n to m), then set the n-th row
% and column to zero.
check = true;
singleDepsCum = [];
while check
    nDependency = sum(abs(HTableRed)>0,2);
    singleDeps = find(nDependency==1);% these are intrinsically increasing
    singleDepsLoc = setdiff(singleDeps,singleDepsCum);
    singleDepsCum = union(singleDeps,singleDepsCum);
    for ind = numel(singleDepsLoc):-1:1
        jnd = find(abs(HTableRed(singleDepsLoc(ind),:))>0);
        if ~isempty(jnd)
            signTable(jnd,singleDepsLoc(ind)) = sign(HTableRed(singleDepsLoc(ind),jnd));
        HTableRed(:,jnd) = HTableRed(:,jnd) + HTableRed(:,singleDepsLoc(ind))*sign(HTableRed(singleDepsLoc(ind),jnd));
        HTableRed(:,singleDepsLoc(ind)) = 0;
        HTableRed(singleDeps(ind),:) = 0;
        end
    end
    check = ~isempty(singleDeps);
end

% Transpose
singleDepsCum = [];
if sum(nDependency)>0
    check = true;
    while check
        nDependency = sum(abs(HTableRed)>0,1);
        singleDeps = find(nDependency==1);
        singleDepsLoc = setdiff(singleDeps,singleDepsCum);
        singleDepsCum = union(singleDeps,singleDepsCum);
        for ind = numel(singleDepsLoc):-1:1
            jnd = find(abs(HTableRed(:,singleDepsLoc(ind)))>0);
            if ~isempty(jnd)
                signTable(jnd,singleDepsLoc(ind)) = sign(HTableRed(jnd,singleDepsLoc(ind)));
                HTableRed(jnd,:) = HTableRed(jnd,:) + HTableRed(singleDepsLoc(ind),:)*sign(HTableRed(jnd,singleDepsLoc(ind)));
                HTableRed(:,singleDepsLoc(ind)) = 0;
                HTableRed(singleDepsLoc(ind),:) = 0;
            end
        end
        check = ~isempty(singleDeps);
    end
end

%play back sign table
inds = find(sum(abs(signTable),1)==0);
sstot = [];
for ind = inds
    ss = zeros(1,size(HTable,2));
    recursiveSignRetrieval(ind,1);
    sstot = cat(1,sstot,ss);
end

% 
% % First, identify rows n with a single dependence on m. Add the n-th column
% % to the m-th column (moving dependencies on n to m), then set the n-th row
% % and column to zero.
% check = true;
% singleDepsCumRow = [];
% singleDepsCumCol = [];
% while check
%     nDependencyRow = sum(abs(HTableRed)>0,2);
%     nDependencyCol = sum(abs(HTableRed)>0,1);
%     singleDepsRow = find(nDependencyRow==1);% these are intrinsically increasing
%     singleDepsCol = find(nDependencyCol==1);% these are intrinsically increasing
%     singleDepsLocRow = setdiff(singleDepsRow,singleDepsCumRow);
%     singleDepsLocCol = setdiff(singleDepsCol,singleDepsCumCol);
%  
%     [valRowCol,colInd] = max(abs(HTableRed(singleDepsLocRow,:))>0,[],2);
%     [valColRow,rowInd] = max(abs(HTableRed(:,singleDepsLocCol))>0,[],1);
%     % find mutual pairs
%     [~,ai,bi] = intersect(singleDepsLocRow,rowInd);
%     ii = find(colInd(ai)==singleDepsLocCol(bi));
%     if ~isempty(ii)
%         for ind = 1:ii
%             singleDepsLocRow(ai(ind))
%             colInd(ai(ind))
%             singleDepsLocCol(bi(ind))
%             rowInd(bi(ind))
%         end
%     end
%     
%     singleDepsCum = union(singleDeps,singleDepsCum);
%     for ind = numel(singleDepsLoc):-1:1
%         jnd = find(abs(HTableRed(singleDepsLoc(ind),:))>0);
%         signTable(jnd,singleDepsLoc(ind)) = sign(HTableRed(singleDepsLoc(ind),jnd));
%         HTableRed(:,jnd) = HTableRed(:,jnd) + HTableRed(:,singleDepsLoc(ind))*sign(HTableRed(singleDepsLoc(ind),jnd));
%         HTableRed(:,singleDepsLoc(ind)) = 0;
%         HTableRed(singleDeps(ind),:) = 0;
%     end
%     check = ~isempty(singleDeps);
% end
% 
% % Transpose
% singleDepsCum = [];
% if sum(nDependency)>0
%     check = true;
%     while check
%         nDependency = sum(abs(HTableRed)>0,1);
%         singleDeps = find(nDependency==1);
%         singleDepsLoc = setdiff(singleDeps,singleDepsCum);
%         singleDepsCum = union(singleDeps,singleDepsCum);
%         for ind = numel(singleDepsLoc):-1:1
%             jnd = find(abs(HTableRed(:,singleDepsLoc(ind)))>0);
%             signTable(jnd,singleDepsLoc(ind)) = sign(HTableRed(jnd,singleDepsLoc(ind)));
%             HTableRed(jnd,:) = HTableRed(jnd,:) + HTableRed(singleDepsLoc(ind),:)*sign(HTableRed(jnd,singleDepsLoc(ind)));
%             HTableRed(:,singleDepsLoc(ind)) = 0;
%             HTableRed(singleDepsLoc(ind),:) = 0;
%         end
%         check = ~isempty(singleDeps);
%     end
% end
% 
% %play back sign table
% inds = find(sum(abs(signTable),1)==0);
% sstot = [];
% for ind = inds
%     ss = zeros(1,size(HTable,2));
%     recursiveSignRetrieval(ind,1);
%     sstot = cat(1,sstot,ss);
% end
% 

%%
function recursiveSignRetrieval(ind,ssign)
    knds = find(abs(signTable(ind,:))>0);
    ss(ind) = ssign;
    if ~isempty(knds)
        for knd = 1:numel(knds)
            recursiveSignRetrieval(knds(knd),signTable(ind,knds(knd))*ssign);
        end
    end
end % function recursiveSignRetrieval()       

end % function findSignDependency

