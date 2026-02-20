% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [TensMat] = Tens2(TrCh,TA)
%TENS2 This is a tensor identical to Tnes1 function, except it uses the
%*squared* surface area of the triangle as wieghts, and finally divedes the
%matrice by the sum of squared areas.
%   As in Tens1 it was designed for unit vectors, but here is can take
%   either these or the coordinates of triangles' centroids
%
% Inputs:
%   TrCh
%   TA
% Outputs:
%   TensMat
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.
dim = size(TrCh,2);

switch dim
    case 3
        Sxx = sum(TA.*(TrCh(:,1).^2));
        Syy = sum(TA.*(TrCh(:,2).^2));
        Szz = sum(TA.*(TrCh(:,3).^2));
        Sxy = sum(TA.*(TrCh(:,1).*TrCh(:,2)));
        Sxz = sum(TA.*(TrCh(:,1).*TrCh(:,3)));
        Syz = sum(TA.*(TrCh(:,2).*TrCh(:,3)));

        TensMat = [Sxx, Sxy, Sxz;...
                   Sxy, Syy, Syz;...
                   Sxz, Syz, Szz;]/sum(TA);
    case 2
        Sxx = sum(TrCh(:,1).^2);
        Syy = sum(TrCh(:,2).^2);
        Sxy = sum(TrCh(:,1).*TrCh(:,2));
        
        TensMat = [Sxx, Sxy;...
                   Sxy, Syy];
end

end

