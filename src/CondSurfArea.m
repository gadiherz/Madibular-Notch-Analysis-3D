% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [CondLn, CondSurf, TriCrvMean, TriCrvMeanSD, TriCrvGaus, TriCrvGausSD] = CondSurfArea(Rotv, RotvRam, RotfRam, TrR, CondLn,CondCntr, RKmax, RKmin)
%CondSurfArea Returns the surface area and sum of absolute mean curvature
% of the area bounded by the Condylar line. It also return a modified
% version of the condilar line so that each line section is a mesh edge.
%By Gadi Herzlinger
%
% Inputs:
%   Rotv
%   RotvRam
%   RotfRam
%   TrR
%   CondLn
%   CondCntr
%   RKmax
%   RKmin
% Outputs:
%   CondLn
%   CondSurf
%   TriCrvMean
%   TriCrvMeanSD
%   TriCrvGaus
%   TriCrvGausSD
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%   It finds all the triangles bounded by the line, and sums their area. It
%   also calculates the sum of absolute mean curvature to surface area
%   ratio to find how "bent" is this region.

%Turn CondLine into a list of edges
CondLnEdg = [CondLn, circshift(CondLn,-1)];
CondLnEdg = CondLnEdg(1:end-1,:);
%find where two adjecent points do not form an edge on the mesh 
NonCon = ~isConnected(TrR,CondLnEdg);

%Use graph approach for finding shortest path
RamEdg = edges(TrR);
RamEdgW = vecnorm(Rotv(RotvRam(RamEdg(:,2)),:)-Rotv(RotvRam(RamEdg(:,1)),:),2,2);
G = graph(RamEdg(:,1),RamEdg(:,2),RamEdgW);  

%fix it by adding intermediate points with direct connections (edges)
%and adding them to the line in the right place. Repeat until the entire
%line is connecting only mesh-adjecent points.
LnPosChk = 0;
k=1;
while any(NonCon)
    LnPos = find(NonCon,1);
    Conn = shortestpath(G,CondLnEdg(LnPos,1),CondLnEdg(LnPos,2));
    Conn = Conn(2:end-1)';
    Conn = Conn(~ismember(Conn,CondLn));
    if isempty(Conn)
        if LnPos == LnPosChk
            k=k+1;
        else
            k=1;
            LnPosChk = LnPos;
        end
        CondLn = [CondLn(1:LnPos-k);CondLn(LnPos+1+k:end)];
    else
        CondLn = [CondLn(1:LnPos);Conn;CondLn(LnPos+1:end)];
    end
    CondLnEdg = [CondLn, circshift(CondLn,-1)];
    CondLnEdg = CondLnEdg(1:end-1,:);
    NonCon = ~isConnected(TrR,CondLnEdg); 
end


%Once the line is alligned with the grid, find the triangles on either side
%of each edge (line segment)
TriBound = cell2mat(edgeAttachments(TrR,CondLnEdg));

%Apply a region growing procedure from the centre point (used to find the
%line in the first place) with a stopping rule using the boundary
%triangles.
StP = nearestNeighbor(TrR,CondCntr);
InCond = cell2mat(vertexAttachments(TrR,StP))';
Add=0;
while ~isempty(Add)
    Add = neighbors(TrR,InCond);
    Add = unique(Add);
    Add = setdiff(Add,InCond);
    Add = setdiff(Add,TriBound);
    InCond = [InCond;Add];
end

TriP1 = Rotv(RotvRam(RotfRam(InCond,1)),:);
TriP2 = Rotv(RotvRam(RotfRam(InCond,2)),:);
TriP3 = Rotv(RotvRam(RotfRam(InCond,3)),:);

TriV1 = TriP2 - TriP1;
TriV2 = TriP3 - TriP2;

TriArea = vecnorm(cross(TriV1,TriV2,2),2,2)./2;
CondSurf = sum(TriArea);
TriArea = [TriArea,TriArea./CondSurf];

TRiPInd = RotfRam(InCond,:);
TriCrvMX = mean(RKmax(TRiPInd),2);
TriCrvMN = mean(RKmin(TRiPInd),2);
TriCrvMean = sum(((TriCrvMX+TriCrvMN)./2).*TriArea(:,2));
TriCrvMeanSD = std(((TriCrvMX+TriCrvMN)./2).*TriArea(:,2));
TriCrvGaus = sum(TriCrvMX.*TriCrvMN.*TriArea(:,2));
TriCrvGausSD = std(TriCrvMX.*TriCrvMN.*TriArea(:,2));
%MeanCurv = sum(abs((RKmax(TRiPInd)+RKmin(TRiPInd))./2),'all');


end

