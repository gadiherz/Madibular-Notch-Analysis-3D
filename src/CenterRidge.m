% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Vatt] = CenterRidge(Tr, RdgEdgs, KDmin)
%CENTERRIDGE Summary of this function goes here
%   Detailed explanation goes here
%
% Inputs:
%   Tr
%   RdgEdgs
%   KDmin
% Outputs:
%   Vatt
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%This section is used to find for each connected ridge point the
%intersection points of adjecent vectors (i.e. edges) and the plane defined
%by the direction of minimum curvature (i.e. the ridge) as a normal to the
%plane. 
%prameter d of line function in vector notation d = (p0-l0).n / l.n Where
%p0 is a point on the plane (the ridge point) l0 is a point on the line
%(the first entry in RdgEdgs) n is the normal vector to the plane (Cmin of
%the ridge point) and l is the direction of the line (the differnce between
%the two ridge points connected by the relevant edge). This is calculated
%in cell array using cellfun to avoid loops.

%***These will then be averaged and apparently snapped to the
%closest ridge point (or vertex on mesh?? assuming the average of these
%intersection points reflects the "middle" of the ridge...)


fRam = Tr.ConnectivityList;
vRam = Tr.Points;

%Find close (4 ring) connected rdgpts and calculate verctor-plane
%intersection only for their vectors (for each pt)
ConRdgPt = unique(RdgEdgs);
clear Vatt

%Finds neighborig ridge points within 6 rings from each ridge points
Vatt = vertexAttachments(Tr,ConRdgPt);

for i =1:2
    Vatt = cellfun(@(x) unique(fRam(x,:)), Vatt, 'uniformoutput', false);
    FlipCells = find(cellfun(@(x) size(x,2), Vatt)~=1);
    for j =1:size(FlipCells,1)
        Vatt{FlipCells(j),1} = Vatt{FlipCells(j),1}';
    end
    Vatt = cellfun(@cell2mat, cellfun(@(x) vertexAttachments(Tr,x)', Vatt, 'UniformOutput', false), 'UniformOutput', false);
end
clear FlipCells
Vatt = cellfun(@(x) unique(fRam(x,:)), Vatt, 'uniformoutput', false);
Vatt = [num2cell(ConRdgPt), Vatt];
Vatt(:,2) = cellfun(@(x) x(ismember(x,ConRdgPt)), Vatt(:,2), 'UniformOutput',false);

%gets indices of vectors passing through these points (in any direction)
Vatt(:,3) = cellfun(@(x) find(any(ismember(RdgEdgs,x),2)), Vatt(:,2), 'UniformOutput', false);
%gets the actual vectors
Vatt(:,4) = cellfun(@(x) vRam(RdgEdgs(x,1),:) - vRam(RdgEdgs(x,2),:),Vatt(:,3),'UniformOutput',false);
%gets a point on the vector (first column of RdgEdgs)
Vatt(:,5) = cellfun(@(x) vRam(RdgEdgs(x,1),:), Vatt(:,3), 'UniformOutput',false);
%gets point on plane (The ridge point's index i.e. Vatt(:,1)), repeats it 
%according to the number of vectors and subtracts the point on line from it
Vatt(:,6) = cellfun(@(x) vRam(x,:), Vatt(:,1), 'UniformOutput',false);
Vatt(:,6) = cellfun(@(x,y) repmat(x,[size(y,1),1]), Vatt(:,6), Vatt(:,5), 'UniformOutput',false);
Vatt(:,6) = cellfun(@(x,y) minus(x,y), Vatt(:,6), Vatt(:,5), 'UniformOutput',false);
%gets plane normal and repeats it according to the number of vectors
Vatt(:,7) = cellfun(@(x) KDmin(x,:),Vatt(:,1), 'UniformOutput',false);
Vatt(:,7) = cellfun(@(x,y) repmat(x,[size(y,1),1]), Vatt(:,7), Vatt(:,5), 'UniformOutput',false);

%dot product of point difference (:,6) and plane normal (:,7) - Note that
%this is allocated to column 6!!!
Vatt(:,6) = cellfun(@(x,y) dot(x,y,2), Vatt(:,6), Vatt(:,7), 'UniformOutput', false);
%Dot product of line vector (:,4) and that normal to the plane (:,7) - Note
%that this is allocated to column 7!!!
Vatt(:,7) = cellfun(@(x,y) dot(x,y,2), Vatt(:,4), Vatt(:,7), 'UniformOutput', false);
%division of the last two values gives the d parameter for calcuclating
%the intersection. This is allocated to column 6 and Vatt is trimmed
%accordingly
Vatt(:,6)  = cellfun(@(x,y) rdivide(x,y), Vatt(:,6), Vatt(:,7), 'UniformOutput', false);
Vatt = Vatt(:,1:6);

%This finds the point of intersection between each of the edges and the
%plane and creates a list containing: (1) The point index from the original 
%connected ridge points.(2) The indices of edges associated woth this point
%(3) The direction of the ridge (min curve) at that point. (4)The mean
%point of intersection between the edge vetors and the plane.
Vatt(:,7) = cellfun(@(x,y,z) x+y.*z, Vatt(:,5), Vatt(:,4), Vatt(:,6), 'UniformOutput', false);
Vatt = [Vatt(:,1),Vatt(:,3),mat2cell(KDmin(cell2mat(Vatt(:,1)),:),ones(size(Vatt,1),1),3),Vatt(:,7)];
Vatt(:,4) = cellfun(@(x) mean(x,1), Vatt(:,4),'UniformOutput', false);
RemNans = cellfun(@(x) any(isnan(x)) | any(isinf(x)), Vatt(:,4));
Vatt=Vatt(~RemNans,:);
% %Reduce ridge points
% Vatt = RedRdgPts(Vatt, RdgEdgs);

%This snaps the mean point coordinate to the closest vertex on the mesh. It
%also updates the vertex indices in Vatt(:,1). Note that the vectors
%pointing in the direction of the ridge are not updated but retained from
%the original ridge points. The loop is used remove duplicate rows. Where 
%there are such cases the vector pointing in the direction of the ridge is
%averaged from all the multiple instances and normalized again.
SnapVer = nearestNeighbor(Tr,cell2mat(Vatt(:,4)));
Vatt(:,4) = mat2cell(vRam(SnapVer,:),ones(size(SnapVer,1),1),3);
Vatt(:,1) = num2cell(SnapVer);
[~, Frow, Drow] = unique(SnapVer,'stable');
for i=1:size(Frow,1)
    DirVecs = cell2mat(Vatt(:,3));
    rows = find(Drow==i);
    if size(rows,1)>1
        Nvec = mean(DirVecs(rows,:));
        Nvec = Nvec/norm(Nvec);
        Vatt(Frow(i),3) = mat2cell(Nvec,1,3);
    end
end
Vatt = Vatt(Frow,:);
clear Frow Drow DirVecs rows Nvec
end

