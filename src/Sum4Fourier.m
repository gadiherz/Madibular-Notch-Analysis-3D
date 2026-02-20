% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Val] = Sum4Fourier(TDist, Coeffs, CuDist, Wn)
%Sum4Fourier sums the terms in the Fourier Transform function
%   It receives the total normalized distance (=1), the coefficients
%   (either X or Y) and the cumulative distance along the curve. It also
%   recives the smoothing factor Wn. It returns the sum of all terms in the equation.
%   It is operated by FourierCoords.m
%
% Inputs:
%   TDist
%   Coeffs
%   CuDist
%   Wn
% Outputs:
%   Val
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

for i = 1:size(Coeffs,1)
    Eq(i,1) = Wn(i)*(Coeffs(i,1)*sin(2*pi/TDist*(i-1)*CuDist)+Coeffs(i,2)*cos(2*pi/TDist*(i-1)*CuDist));
end

Val = sum(Eq);

end

