% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [RamSymVec, RMat] = RamSymPln(RvN)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%
% Inputs:
%   RvN
% Outputs:
%   RamSymVec
%   RMat
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.
MatNV = [sum(RvN(:,1).^2),sum(RvN(:,1).*RvN(:,2)),sum(RvN(:,1).*RvN(:,3));...
    sum(RvN(:,2).*RvN(:,1)),sum(RvN(:,2).^2),sum(RvN(:,2).*RvN(:,3));...
    sum(RvN(:,3).*RvN(:,1)),sum(RvN(:,3).*RvN(:,2)),sum(RvN(:,3).^2)];
[RMat,~] = eig((MatNV));
RamSymVec = (RMat*[0,0,1]');
end