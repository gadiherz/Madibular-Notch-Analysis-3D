% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Ang] = RamCondAng(RvN,CondRMat)
%RamCondAng Returns the angle between the normal vector to the ramus
%symmetry plane and the corresponding eigenvector of the condyle (the
%"lateral vector")
%
% Inputs:
%   RvN
%   CondRMat
% Outputs:
%   Ang
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

% calculate ramus coordinate system based on the surfaces' normal vectors
MatRam = [sum(RvN(:,1).^2),sum(RvN(:,1).*RvN(:,2)),sum(RvN(:,1).*RvN(:,3));...
    sum(RvN(:,2).*RvN(:,1)),sum(RvN(:,2).^2),sum(RvN(:,2).*RvN(:,3));...
    sum(RvN(:,3).*RvN(:,1)),sum(RvN(:,3).*RvN(:,2)),sum(RvN(:,3).^2)];
[eVecsRam,eValsRam] = eig((MatRam));
%For the ramus the vector is the one with the highest eignvalue
[~,SelVec] = find(abs(eValsRam)==max(abs(eValsRam),[],"all"));
RamVec = eVecsRam(:,SelVec);

%For the condyle it is more difficult to identify the correct vector. In the
%given rotation matrix built in the CondLnPrjC3 function it is either the
%first of third vector (X/Z). Here both are tested and the one with the
%lower angle is the correct one.

Ang(1) = acos(dot(RamVec',CondRMat(:,1)',2));
Ang(2) = acos(dot(RamVec',CondRMat(:,3)',2));
%rearrange in 0<ang<0.5*pi
Ang = abs(abs(Ang-0.5*pi)-0.5*pi); 
Ang = min(Ang);
end