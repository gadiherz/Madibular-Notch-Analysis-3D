% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Rdgnss] = CalcRdgnss(Rotv,Kmax,Kmin,KDmax)
%CALCRDGNSS Returns the "ridgness" value of each vertex given vertex
%coordinate maximum principal curvature and its direction
%   The function finds the neighbors of a vertex within a given range of 
%   the slope direction and checks the relation of the vertex's curvature
%   to that of its neighbors - a vertex whoes curvature is higher than all
%   of its relevant neighbors will have a high ridgness value.
%
% Inputs:
%   Rotv
%   Kmax
%   Kmin
%   KDmax
% Outputs:
%   Rdgnss
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%Get 200 nearest neighbors of each point
NSz = 100;
Vatt = knnsearch(Rotv,Rotv,'k',NSz);

%Get normalized vector between each vector and point (3D vectorization storage)
Vdiff = zeros(size(Vatt,1),size(Vatt,2)-1,3);
Vdiff(:,:,1) = reshape(Rotv(Vatt(:,2:end),1),size(Vatt,1),[]) - Rotv(Vatt(:,1),1);
Vdiff(:,:,2) = reshape(Rotv(Vatt(:,2:end),2),size(Vatt,1),[]) - Rotv(Vatt(:,1),2);
Vdiff(:,:,3) = reshape(Rotv(Vatt(:,2:end),3),size(Vatt,1),[]) - Rotv(Vatt(:,1),3);

%find distance and normalize
Vdist = vecnorm(Vdiff,2,3);
Vdiff = Vdiff./Vdist;

%Find angle between vectors above and direction of max curvature
KDmax3 = repmat(pagetranspose(shiftdim(KDmax,-1)),1,NSz-1,1);
Vang = acos(dot(Vdiff,KDmax3,3));
%mark and count those to include
Vinc1 = Vang<=deg2rad(45); 
Vinc2 = Vang>=deg2rad(135);
Vcnt1 = sum(Vinc1,2);
Vcnt2 = sum(Vinc2,2);
%difference between max and the absolute of min curvature of the point and its neighbors
Kmax3 = Kmax(Vatt(:,2:end));%-Kmin(Vatt(:,2:end));
KmaxP3 = repmat(Kmax,1,NSz-1);%-Kmin !!put back

%get difference relative to distance and normalize to number of neighbors
Vres1 = zeros(size(Vang));
Vres2 = zeros(size(Vang));
%Vres1(Vinc1) = (KmaxP3(Vinc1)-0.95*Kmax3(Vinc1))./(Vdist(Vinc1));%+1).^2;
%Vres2(Vinc2) = (KmaxP3(Vinc2)-0.95*Kmax3(Vinc2))./(Vdist(Vinc2));%+1).^2;

Vres1(Vinc1) = (KmaxP3(Vinc1)>0.95*Kmax3(Vinc1))./(Vdist(Vinc1));%+1).^2;
Vres2(Vinc2) = (KmaxP3(Vinc2)>0.95*Kmax3(Vinc2))./(Vdist(Vinc2));%+1).^2;


Sd1 = sum(Vres1,2)./Vcnt1;
Sd2 = sum(Vres2,2)./Vcnt2;
Sd1(isnan(Sd1)) = 0;
Sd2(isnan(Sd2)) = 0;
Rdgnss = (Sd1 + Sd2 - abs(Sd1 - Sd2));
Rdgnss(Rdgnss>0) = Kmax(Rdgnss>0);
%Rdgnss = Sd1+Sd2;

%clean outliers and return
%Rdgnss = filloutliers(Rdgnss,'clip');
end

