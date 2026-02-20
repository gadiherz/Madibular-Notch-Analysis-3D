% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [EdgePts,TanAng] = FindCondEdgPts(RefPt,PrjRotv,RotvRam,InV,RvN,RotfCnd)
%FindCondEdgPts Returns the edge points on a condyle
%   The function gets a reference points which is "in the middle" (on the
%   plane) and "above" the condyle as well as the rotated vertex
%   coordinates, ramus and condyle (based on notch minimum) indices,
%   normals and a approproatly indexed triangle list. It returns the
%   indices (refering to RotvRam) of edge points whose normals are ~1/4*pi
%   to the vectors between them and the reference point (tangential).
%
% Inputs:
%   RefPt
%   PrjRotv
%   RotvRam
%   InV
%   RvN
%   RotfCnd
% Outputs:
%   EdgePts
%   TanAng
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%vectors, directions and angles between all vertices on condyle and ref points
TanVec = PrjRotv(RotvRam(InV),:)-RefPt;
TanVecDir = TanVec./vecnorm(TanVec,2,2);
TanAng = acos(dot(RvN(InV,:),TanVecDir,2));
%Addition
TanAngCrd = acos(dot(TanVecDir,repmat([0,-1,0],size(TanVecDir,1),1),2));

%Keep only edgepoints and recalc all, keep magnitude
EdgePts = InV(TanAng>deg2rad(85) & TanAng<deg2rad(95));

%Addition - changed the condition of edge point to relay on the coordinate
%vector and not justthe normal
EdgePts = InV(TanAng>deg2rad(85) & TanAng<deg2rad(95));%...
%    & TanAngCrd>deg2rad(40) & TanAngCrd<deg2rad(50)); 
EdgePtsInd = find(TanAng>deg2rad(85) & TanAng<deg2rad(95));%...
%    & TanAngCrd>deg2rad(85) & TanAngCrd<deg2rad(95));
CondCrds = RvN(EdgePts,:);
TanVec = PrjRotv(RotvRam(EdgePts),:)-RefPt;
TanVecMag = vecnorm(TanVec,2,2);
TanVecDir = TanVec./TanVecMag;
% TanAng = acos(dot(CondCrds,TanVecDir,2));

%Make sure that the vectors between edge and the reference points do not
%intersect any triangle on the condyle

%Build rotation matrix for each edgepoint 
Rvec = cross(TanVecDir,repmat([0,1,0],size(TanVecDir,1),1),2);
Rang = acos(dot(TanVecDir,repmat([0,1,0],size(TanVecDir,1),1),2));
Rmat = zeros(3,3,size(EdgePts,1));
for i=1:size(EdgePts,1)
    Rmat(:,:,i) = RMatAxAng(Rvec(i,:),-Rang(i));
end
% %Get vertex coords
% Vs = PrjRotv(RotvRam(InV),:)-RefPt;
% 
% %Rotate so that edgepoint is origin and vector between it and reference
% %point is Y positive
% NonInd=1;
% for i=1:size(EdgePtsInd,1)
%     VsRot = (Rmat(:,:,i)*Vs')';
%     %Check if any triangle *intersects* Y, remove triangles touching Y
%     OpTri = find(~(all([VsRot(RotfCnd(:,1),1),VsRot(RotfCnd(:,2),1),VsRot(RotfCnd(:,3),1)]<0,2) | all([VsRot(RotfCnd(:,1),1),VsRot(RotfCnd(:,2),1),VsRot(RotfCnd(:,3),1)]>0,2))...
%         & ~(all([VsRot(RotfCnd(:,1),3),VsRot(RotfCnd(:,2),3),VsRot(RotfCnd(:,3),3)]<0,2) | all([VsRot(RotfCnd(:,1),3),VsRot(RotfCnd(:,2),3),VsRot(RotfCnd(:,3),3)]>0,2)));
%     OpTri = setdiff(OpTri,OpTri(any(RotfCnd(OpTri,:)==EdgePtsInd(i),2)));
%     if ~isempty(OpTri)
%         RelPos=[];
%         for j=1:numel(OpTri)
%             RelPos(j) = mean(VsRot(RotfCnd(OpTri(j),:),2));
%         end
%         %Remove triangle whose Y in negative (below edge point) or above
%         %magnitude
%         OpTri = OpTri(RelPos>0 & RelPos<TanVecMag(i));
%         IntAng = [];
%         for j=1:size(OpTri,1)
%             MVec = mean(RvN(InV(RotfCnd(OpTri(j),:)),:))./vecnorm(mean(RvN(InV(RotfCnd(OpTri(j),:)),:)),2,2);
%             IntAng(j) = acos(dot(MVec,TanVecDir(i,:),2));
%         end
%         OpTri = OpTri(IntAng < deg2rad(80) | IntAng > deg2rad(100));
%     end
%     %If it intersects a triangle remove from edgepoints list
%     if ~isempty(OpTri)
%         NonEdg(NonInd) = EdgePts(i);
%         NonInd = NonInd+1;
%     end
% end
% EdgePts = setdiff(EdgePts,NonEdg);
end