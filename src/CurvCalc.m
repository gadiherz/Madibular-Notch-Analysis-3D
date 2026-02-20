% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [K] = CurvCalc(Dists,Wn,CoeffsX, CoeffsY)
%CurvCalc Calculated the curvature at each point along the curve as a
%function of the cumulative distance along the curve
%   It recives the cumulative distance along the curve 
%   the X and Y coefficients and the smoothing weight
%   parameter vector
%
% Inputs:
%   Dists
%   Wn
%   CoeffsX
%   CoeffsY
% Outputs:
%   K
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%for each distance in the cumulative distances vector it calculates the the
%summed (S) terms of the first (dxds; dyds) and second (d2xds2; d2yds2) 
%derivatives of the original functions. 
TDist = Dists(end);
for j = 1:size(Dists,1)
    Sdxds=0;
    Sdyds=0;
    Sd2xds2=0;
    Sd2yds2=0;
    for i = 1:size(CoeffsX,1)
        dxds = -2*pi*(i-1)*Wn(i)*(CoeffsX(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsX(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;
        dyds = -2*pi*(i-1)*Wn(i)*(CoeffsY(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsY(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;
        d2xds2 = -4*pi.^2*(i-1).^2*Wn(i)*(CoeffsX(i,1)*sin(2*pi*(i-1)*(Dists(j)./TDist)) + CoeffsX(i,2)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist.^2;
        d2yds2 = -4*pi.^2*(i-1).^2*Wn(i)*(CoeffsY(i,1)*sin(2*pi*(i-1)*(Dists(j)./TDist)) + CoeffsY(i,2)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist.^2;
%         dxds = -2*pi*i*Wn(i)*(CoeffsX(i,2)*sin(2*pi*i*Dists(j)) - CoeffsX(i,1)*cos(2*pi*i*Dists(j)))./;
%         dyds = -2*pi*i*Wn(i)*(CoeffsY(i,2)*sin(2*pi*i*Dists(j)) - CoeffsY(i,1)*cos(2*pi*i*Dists(j)));
%         d2xds2 = -4*pi.^2*i.^2*Wn(i)*(CoeffsX(i,1)*sin(2*pi*i*Dists(j)) + CoeffsX(i,2)*cos(2*pi*i*Dists(j)));
%         d2yds2 = -4*pi.^2*i.^2*Wn(i)*(CoeffsY(i,1)*sin(2*pi*i*Dists(j)) + CoeffsY(i,2)*cos(2*pi*i*Dists(j)));
        Sdxds = Sdxds+dxds;
        Sdyds = Sdyds+dyds;
        Sd2xds2 = Sd2xds2+d2xds2;
        Sd2yds2 = Sd2yds2+d2yds2;
        %Kt = Kt + (2*(dxds-dyds)*d2xds2*d2yds2)/(dxds.^2+dyds.^2).^2;
    end
% The curvature at each point is calculated as
% (y'(c)*x''(c)-x'(c)*y''(c))/(x'(c)^2+y'(c)^2)
% This is in fact the second derivative of the function XY fourier
% transform as the first is the angle of the tangent at each point (theta)
% calculated by TanCalc.m
K(j,1) = -(Sdyds*Sd2xds2-Sdxds*Sd2yds2)/(Sdxds.^2+Sdyds.^2);
end
if mean(K)<0
    K=-K;
end
end

