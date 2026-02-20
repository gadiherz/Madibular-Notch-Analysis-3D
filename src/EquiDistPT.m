% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [B2, CuDist] = EquiDistPT(Bcoords, OpenCurv)
%
% Inputs:
%   Bcoords
%   OpenCurv
% Outputs:
%   B2
%   CuDist
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

% EquiDistPT returns equdistantly spaced points and their cumulative distance
% along the curve sampled by Bcoords. It is used by FourierCoefs for the
% following fitting of the Fourier series. It makes sure that the number of
% points is even. 

%Find the curve length
PDist = vecnorm(Bcoords(2:end,:)-Bcoords(1:end-1,:),2,2);
TDist = sum(PDist);

%The required distance is the mean point to point distance as long as it is
%higher than a minimum threshold
RqDist = TDist/(size(Bcoords,1)-1);

if RqDist < 0.75
    RqDist = 0.75;
end

%make sure the number of points will be even (note - for closed curves this
%includes the firts point twice (first&last). Gets the new point to point
%distance
RqPts = fix(TDist./RqDist);
if rem(RqPts,2) == 0 
    RqPts = RqPts-1;
end
RqDist = TDist/RqPts;

%Point placing loop. The points are placed at equal distance along the
%ORIGINAL curve (B1). Note that the Euclidean distance between the sampled
%points is not uniform as thos point sample the original curve. 
B1=Bcoords;
B2 = B1(1,:);
n=2;

while size(B2,1) < RqPts
    DistTrk = vecnorm(B1(n,:) - B2(end,:),2,2);
    if sum(DistTrk) < RqDist
        while sum(DistTrk) < RqDist && n < size(B1,1)
            n=n+1;
            DistTrk = [DistTrk; vecnorm(B1(n-1,:) - B1(n,:),2,2)];
        end
        DirVec = (B1(n,:) - B1(n-1,:))./vecnorm(B1(n,:) - B1(n-1,:),2,2);
        DistVec = RqDist-sum(DistTrk(1:end-1));
        B2(end+1,:) = B1(n-1,:) + DirVec*DistVec;
    else
        DirVec = (B1(n,:) - B2(end,:))./vecnorm(B1(n,:) - B2(end,:),2,2);
        B2(end+1,:) = B2(end,:) + DirVec*RqDist;
    end
end
%Closes the curve
if OpenCurv
    B2 = [B2;B1(end,:)];
else
    B2 = [B2;B2(1,:)];
end

%The cunulative distance 
CuDist = [0; cumsum(repmat(RqDist,size(B2,1)-1,1))];


%Graphics for testing
    
% Fig1=figure();
% Ax1 = axes(Fig1);
% hold(Ax1,'on');
% BcSc = scatter3(Bcoords(:,1),Bcoords(:,2),Bcoords(:,3),'ob');
% BcPl = plot3(Bcoords(:,1),Bcoords(:,2),Bcoords(:,3),'b');
% B2Sc = scatter3(B2(:,1),B2(:,2),B2(:,3),'ok','filled');
% B2Pl = plot3(B2(:,1),B2(:,2),B2(:,3),'k');
end