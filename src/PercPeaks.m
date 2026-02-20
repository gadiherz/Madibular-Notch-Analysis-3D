% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Peaks] = PercPeaks(Peaks,K,NewDist)
%PercPeaks optimizes the peaks locations after increasing the resolution of K
%   The original input Peaks var is calculated on a low resolution curve
%   (Neff=10). The input K and NewDist are calculated after a second iteration
%   in higher res (Neff=45). This function updates the index, normalized
%   position and magnitude of the peak on the higher res curve by seraching the
%   corresponding max/min K in a range of +-3.5% form the original peak. if the
%   peak doesn't correspond to the sign of curvature the values of that
%   peak remains the same.
%
% Inputs:
%   Peaks
%   K
%   NewDist
% Outputs:
%   Peaks
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.
TDist = NewDist(end);
NewDist = NewDist./TDist;
Peaks(:,2) = Peaks(:,2)./TDist;
for i =1:size(Peaks,1)
    if Peaks(i,2) - 0.015 < 0 
        TmpRangeInd = find((NewDist(:,1) > 0 & NewDist(:,1) < Peaks(i,2) + 0.015)...
            | (NewDist(:,1) > Peaks(i,2) - 0.015+1));
    elseif Peaks(i,2) + 0.015 > 1 
        TmpRangeInd = find((NewDist(:,1) > Peaks(i,2) - 0.015 & NewDist(:,1) < 1)...
            | (NewDist(:,1) < Peaks(i,2) + 0.015-1));
    else
        TmpRangeInd = find((NewDist(:,1) > Peaks(i,2) - 0.015 & NewDist(:,1) < Peaks(i,2) + 0.015));
    end
    if Peaks(i,3) < 0 
        TmpPeaksN(1,1) = TmpRangeInd(find(K(TmpRangeInd,1) == min(K(TmpRangeInd,1))));
        TmpPeaksN(1,2) = NewDist(TmpPeaksN(1,1)).*TDist;
        TmpPeaksN(1,3) = K(TmpPeaksN(1,1));
        if TmpPeaksN(1,3) < 0
            OK = true;
        else
            OK = false;
        end
    else
        TmpPeaksN(1,1) = TmpRangeInd(find(K(TmpRangeInd,1) == max(K(TmpRangeInd,1))));
        TmpPeaksN(1,2) = NewDist(TmpPeaksN(1,1)).*TDist;
        TmpPeaksN(1,3) = K(TmpPeaksN(1,1));
        if TmpPeaksN(1,3) > 0
            OK = true;
        else
            OK = false;
        end
    end
    if OK
        Peaks(i,:) = TmpPeaksN(1,:);
    end
    clear TmpRangeInd TmpPeaksN OK
end

