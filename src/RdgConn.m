% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [RdgEdgs] = RdgConn(Tr,RdgPt,RdgDir)
%RDGCONN Summary of this function goes here
%   Detailed explanation goes hereCo65
%
% Inputs:
%   Tr
%   RdgPt
%   RdgDir
% Outputs:
%   RdgEdgs
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

f = Tr.ConnectivityList;
v = Tr.Points;

Vatt = vertexAttachments(Tr,RdgPt);

for i =1:4
    Vatt = cellfun(@(x) unique(f(x,:)), Vatt, 'uniformoutput', false);
    FlipCells = find(cellfun(@(x) size(x,2), Vatt)~=1);
    for j =1:size(FlipCells,1)
        Vatt{FlipCells(j),1} = Vatt{FlipCells(j),1}';
    end
    Vatt = cellfun(@cell2mat, cellfun(@(x) vertexAttachments(Tr,x)', Vatt, 'UniformOutput', false), 'UniformOutput', false);
end
Vatt = cellfun(@(x) unique(f(x,:)), Vatt, 'uniformoutput', false);
clear FlipCells

%Find which of these neighbors is also a ridge point
Vatt = cellfun(@setdiff, Vatt, num2cell(RdgPt), 'UniformOutput', false);
Vatt(1:end,2) = {RdgPt};
Vatt(:,1) = cellfun(@(x,y) x(find(ismember(x,y))), Vatt(:,1),Vatt(:,2), 'UniformOutput', false);
Vatt = Vatt(:,1);

%finds the vector between each ridge point and all its neighboring ridge
%points
Vatt(:,2) = cellfun(@(x,y) minus(v(x,:),v(y,:)), Vatt(:,1), num2cell(RdgPt) ,'UniformOutput', false);

%Calls the ridge direction vector (RdgDir) of RdgPt(i) and repeats it according to
%the number of neighbors (required by the following)
Vatt(:,3) = cellfun(@(x) RdgDir(x,:),num2cell(RdgPt),'UniformOutput',false);
Vatt(:,4) = cellfun(@(x) size(x,1),Vatt(:,1),'UniformOutput',false);
Vatt(:,5) = {1};
Vatt(:,3) = cellfun(@repmat, Vatt(:,3), Vatt(:,4),Vatt(:,5),'UniformOutput',false);
Vatt = Vatt(:,1:3);

%Finds the angle between the ridge direction vector of the examined point
%RdgPt(i) and those of its neighbors (arc cosine of the dot product between 
%them, devided by the product of their lengths 
Vatt(:,4) = cellfun(@dot, Vatt(:,2),Vatt(:,3),num2cell(repmat(2,size(Vatt,1),1)),'UniformOutput',false);
Vatt(:,5) = cellfun(@sqrt, cellfun(@dot, Vatt(:,2),Vatt(:,2),num2cell(repmat(2,size(Vatt,1),1)),'UniformOutput',false),'UniformOutput',false);
Vatt(:,2) = cellfun(@rdivide, Vatt(:,4), Vatt(:,5),'UniformOutput',false);
Vatt(:,2) = cellfun(@acos, Vatt(:,2),'UniformOutput',false);
Vatt = Vatt(:,1:2);

%Leaves only neighboring points that fall within 5 degrees from the
%RdgDir (direction of the ridge), than finds the
%two closest ones on each direction (min - close to 0 and max close to pi)
Vatt(:,3) = cellfun(@(x) find(x<=deg2rad(10) | x>=deg2rad(170)), Vatt(:,2),'UniformOutput',false);
Vatt(:,1) = cellfun(@(x,y) x(y), Vatt(:,1), Vatt(:,3),'UniformOutput',false);
Vatt(:,2) = cellfun(@(x,y) x(y), Vatt(:,2), Vatt(:,3),'UniformOutput',false);
Vatt = Vatt(:,1:2);

%Vatt(:,3) = cellfun(@(x) find(x == max(x) | x == min(x)), Vatt(:,2), 'UniformOutput', false);
Vatt(:,3) = cellfun(@(x) find((x-0.5*pi == max(x-0.5*pi) & x-0.5*pi > 0) |  (x-0.5*pi == min(x-0.5*pi) & x-0.5*pi < 0)),...
  Vatt(:,2), 'UniformOutput', false);

%New dataset including the ridge points indices that have two such
%neighbors, and the number of their neighbors
Vatt = [num2cell(RdgPt(cellfun(@(x) size(x,1)>=1, Vatt(:,3)))), Vatt(cellfun(@(x) size(x,1)>=1, Vatt(:,3)),1), Vatt(cellfun(@(x) size(x,1)>=1, Vatt(:,3)),3)];
Vatt(:,2) = cellfun(@(x,y) x(y), Vatt(:,2), Vatt(:,3),'UniformOutput', false);
Vatt = Vatt(:,1:2);

%Creates a list of ridge connection between two relevant ridge points. Note
%that each ridge point may have more than two connections so the ridge it
%marked by multiple "subparallel" lines
Vatt(cellfun(@(x) size(x,1)==2, Vatt(:,2)),1) = cellfun(@(x) [x;x], Vatt(cellfun(@(x) size(x,1)==2, Vatt(:,2)),1),'UniformOutput',false);
RdgEdgs = cell2mat(Vatt);


end

