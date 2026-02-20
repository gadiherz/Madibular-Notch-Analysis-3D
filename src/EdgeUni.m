% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [NewRdgPts,NewRdgEdgs] = EdgeUni(Tr, RdgEdgs, KDmin)
%
% Inputs:
%   Tr
%   RdgEdgs
%   KDmin
% Outputs:
%   NewRdgPts
%   NewRdgEdgs
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

fRam = Tr.ConnectivityList;
vRam = Tr.Points;

RdgEdgs = sort(RdgEdgs,2);
RdgEdgs = unique(RdgEdgs,"rows");
ConRdgPt = unique(RdgEdgs);

NewRdgPts = nan(size(ConRdgPt));
NewRdgEdgs = nan(size(RdgEdgs));
Tested = nan(size(ConRdgPt));
while ~isempty(ConRdgPt)
    NxtPt = find(~isnan(NewRdgPts),1,'last')+1;
    [NxtRdg,~] = find(~isnan(NewRdgEdgs),1,'last');
    NxtRdg = NxtRdg+1;
    NxtTst = find(~isnan(Tested),1,'last')+1;
    
    CurPt = ConRdgPt(1);
    %Calculate d parameter of the plane defined by the current point and
    %the ridge direction (i.e. direction of min curvature)
    d = sum(vRam(CurPt,:).*KDmin(CurPt,:),2);
    
    %trim ridge point list
    ConRdgPt = ConRdgPt(2:end);
    
    %Finds neighborig ridge points within 6 rings from each ridge points
    Vatt = vertexAttachments(Tr,CurPt);
    
    for i =1:6
        Vatt = cellfun(@(x) unique(fRam(x,:)), Vatt, 'uniformoutput', false);
        FlipCells = find(cellfun(@(x) size(x,2), Vatt)~=1);
        for j =1:size(FlipCells,1)
            Vatt{FlipCells(j),1} = Vatt{FlipCells(j),1}';
        end
        Vatt = cellfun(@cell2mat, cellfun(@(x) vertexAttachments(Tr,x)', Vatt, 'UniformOutput', false), 'UniformOutput', false);
    end
    clear FlipCells
    Vatt = cell2mat(cellfun(@(x) unique(fRam(x,:)), Vatt, 'uniformoutput', false));
    
    %Neighboring point indices which are also ridge points
    Vatt = Vatt(ismember(Vatt,ConRdgPt));
    
    %Get all edges associated with neighboring ridge points and their
    %endpoints
    NeiEdgsInd = find(any(ismember(RdgEdgs,[CurPt;Vatt]),2));
    NeiEdgs = RdgEdgs(NeiEdgsInd,:);
    EndPts = [vRam(NeiEdgs(:,1),:),vRam(NeiEdgs(:,2),:)];
    
    %Calculate the position of the endpoints relative to the plane
    EndPtsPos = [sum(EndPts(:,1:3).*repmat(KDmin(CurPt,:),size(EndPts,1),1),2)-d,...
        sum(EndPts(:,4:6).*repmat(KDmin(CurPt,:),size(EndPts,1),1),2)-d];
    %Find which edges intersect the plane or have an endpoint on the plane
    RelEdgs = sign(EndPtsPos(:,1)) ~= sign(EndPtsPos(:,2)) | any(EndPtsPos==0,2);
    
    %Keep only relevant edges
    NeiEdgsInd = NeiEdgsInd(RelEdgs,:);
    NeiEdgs = NeiEdgs(RelEdgs,:);
    EndPts = EndPts(RelEdgs,:);
    EndPtsPos = EndPtsPos(RelEdgs,:);
    
    if numel(NeiEdgsInd)>7
        %Sort Edges and endpoints so that those on each side of the plane are
        %assigned to different columns
        [~,I] = sort(EndPtsPos,2);
        EndPts(I(:,1)==2,:) = [EndPts(I(:,1)==2,4:6),EndPts(I(:,1)==2,1:3)];
        NeiEdgs(I(:,1)==2,:) = [NeiEdgs(I(:,1)==2,2),NeiEdgs(I(:,1)==2,1)];
        RdgEdgs(NeiEdgsInd,:) = NeiEdgs;
        clear I EndPtsPos
        
        %Unite all points on either side of the plane to a single point and
        %snap it to the nearest vertex
        NewCrd = mean(EndPts,1);
        SnapVer = nearestNeighbor(Tr,[NewCrd(1:3);NewCrd(4:6)])';
        
        %change point indices which were averaged and snapped in original edge 
        % list. Note that this can include edges which were not used in the
        % current run.
        RdgEdgs(ismember(RdgEdgs,NeiEdgs(:,1))) = SnapVer(1);
        RdgEdgs(ismember(RdgEdgs,NeiEdgs(:,2))) = SnapVer(2);
        
        %Remove duplicate (rows in RdgEdgs are sorted)
        RdgEdgs = unique(RdgEdgs,'rows','stable');
    
        %add new values to new list
        if isempty(NxtPt)
            NewRdgPts(1:3) = [SnapVer(1),CurPt,SnapVer(2)]; 
            NewRdgEdgs(1:2,:) = [SnapVer(1),CurPt; CurPt,SnapVer(2)];
            Tested(1) = CurPt;
        else
            NewRdgPts(NxtPt:NxtPt+2) = [SnapVer(1),CurPt,SnapVer(2)];
            NewRdgEdgs(NxtRdg:NxtRdg+1,:) = [SnapVer(1),CurPt; CurPt,SnapVer(2)];
            Tested(NxtTst) = CurPt;
        end
    else
        if isempty(NxtTst)
            NxtTst =1;
        end
        Tested(NxtTst) = CurPt;
    end
        
    %take the first ridge pointon the list (random)
    ConRdgPt = unique(RdgEdgs);
    ConRdgPt = setdiff(ConRdgPt,Tested);
end

NewRdgPts = NewRdgPts(~isnan(NewRdgPts));
[r,~] = find(all(isnan(NewRdgEdgs),2),1,'first');
NewRdgEdgs = NewRdgEdgs(1:r-1,:);


end