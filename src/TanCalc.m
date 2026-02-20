% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Theta] = TanCalc(Dists,Wn,CoeffsX, CoeffsY)
%TanCalc Calculated the tangent angle (theta) at each point along the curve as a
%function of the cumulative distance along the curve
%   It recives the cumulative distance along the curve (must be
%   normalized!!!) the X and Y coefficients and the smoothing weight
%   parameter vector
%
% Inputs:
%   Dists
%   Wn
%   CoeffsX
%   CoeffsY
% Outputs:
%   Theta
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%for each distance in the cumulative distances vector it calculates the the
%summed (S) terms of the first (dxds; dyds) derivatives of the original
%function
TDist = Dists(end);
for j = 1:size(Dists,1)
    Sdxds=0;
    Sdyds=0;
    for i = 1:size(CoeffsX,1)
        dxds = -2*pi*(i-1)*Wn(i)*(CoeffsX(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsX(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;
        dyds = -2*pi*(i-1)*Wn(i)*(CoeffsY(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsY(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;

%         dxds = -2*pi*i*Wn(i)*(CoeffsX(i,2)*sin(2*pi*i*Dists(j)) - CoeffsX(i,1)*cos(2*pi*i*Dists(j)));
%         dyds = -2*pi*i*Wn(i)*(CoeffsY(i,2)*sin(2*pi*i*Dists(j)) - CoeffsY(i,1)*cos(2*pi*i*Dists(j)));
        Sdxds = Sdxds+dxds;
        Sdyds = Sdyds+dyds;
    end
%The theta (in radians) is the Arctan of the quontient of the derivatives
%of y and x functions.
Theta(j,1) = atan(Sdyds/Sdxds);
end
end

