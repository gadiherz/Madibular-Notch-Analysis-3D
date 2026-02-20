% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [PIndPos,PIndNeg] = FindPeaks(K, NewDist,Closed)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% Inputs:
%   K
%   NewDist
%   Closed
% Outputs:
%   PIndPos
%   PIndNeg
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

clear Peaks PIndPos PIndNeg RotPC SecInd SecPeak

if Closed
    K=K(1:end-1);
end

%Ceckes the same for points 2 -> end-1. Last point is checked separatly.
for i = 2:size(K,1)-1
    if K(i) > K(i-1) && K(i) > K(i+1) && K(i) > 0
        Peaks(i,1) = true;
        Peaks(i,2) = false;
    elseif K(i) < K(i-1) && K(i) < K(i+1) && K(i) < 0
        Peaks(i,1) = false;
        Peaks(i,2) = true;
    else
        Peaks(i,1) = false;
        Peaks(i,2) = false;
    end
end

if Closed
    %Checkes if first point is a peak - namely of it is higher or lower then
    %both its predecessor (in this case the last point - curve is closed) and
    %its follower. Also a poisitive/negative peak need not only to be the
    %highest/lowest in its neighborhood but also to actually be
    %positive/negative (<> 0)
    if K(1) > K(end) && K(1) > K(2) && K(1) > 0
        Peaks(1,1) = true;
        Peaks(1,2) = false;
    elseif K(1) < K(end) && K(1) < K(2) && K(1) < 0
        Peaks(1,1) = false;
        Peaks(1,2) = true;
    else
        Peaks(1,1) = false;
        Peaks(1,2) = false;
    end


    %check the same for the last point. As the curve is closed its follower is
    %point one again.
    if K(end) > K(end-1) && K(end) > K(1) && K(end) > 0
        Peaks(end+1,1) = true;
        Peaks(end,2) = false;
    elseif K(end) < K(end-1) && K(end) < K(1) && K(end) < 0
        Peaks(end+1,1) = false;
        Peaks(end,2) = true;
    else
        Peaks(end+1,1) = false;
        Peaks(end,2) = false;
    end
else
    Peaks(1,1) = false;
    Peaks(1,2) = false;
    Peaks(end+1,1) = false;
    Peaks(end,2) = false;
end

PIndPos = find(Peaks(:,1));
PIndPos(:,2) = NewDist(PIndPos(:,1));
PIndPos(:,3) = K(PIndPos(:,1));
PIndNeg = find(Peaks(:,2));
PIndNeg(:,2) = NewDist(PIndNeg(:,1));
PIndNeg(:,3) = K(PIndNeg(:,1));

[~, Ind] = sort(PIndPos(:,3),'descend');
PIndPos = PIndPos(Ind,:);
clear Ind
[~, Ind] = sort(PIndNeg(:,3));
PIndNeg = PIndNeg(Ind,:);
clear Ind

end

