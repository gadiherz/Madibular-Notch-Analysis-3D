% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function RMat = RMatAxAng(Ax,ang)
%RMatAxAng Returns a 3D rotation matrix given an axis and an angle
%   Detailed explanation goes here
%
% Inputs:
%   Ax
%   ang
% Outputs:
%   RMat
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

if vecnorm(Ax,2,2) ~= 1
    Ax = Ax./vecnorm(Ax,2,2);
end
s = cos(-ang);
t = sin(-ang);

RMat = [s+Ax(1).^2.*(1-s), Ax(1).*Ax(2).*(1-s)-Ax(3).*t, Ax(1).*Ax(3).*(1-s)+Ax(2).*t;...
    Ax(2).*Ax(1).*(1-s)+Ax(3).*t, s + Ax(2).^2.*(1-s), Ax(2).*Ax(3).*(1-s)-Ax(1).*t;...
    Ax(3).*Ax(1).*(1-s)-Ax(2).*t, Ax(3).*Ax(2).*(1-s)+Ax(1).*t, s+Ax(3).^2.*(1-s)];


end