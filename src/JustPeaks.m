% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Peaks] = JustPeaks(K, NewDist, Closed)
%JustPeaks finds the curvature peaks (K)
%   The output if a var containing the index of the peak (:,1), the
%   normalized length along the curve (:,2), and its magnitude (:,3). The
%   output is given only if it can reach a configuration in which rows 1 &
%   3 are both negative, while all other are positive (the sock paradigm),
%   otherwise the Peaks variable is returned empty
%
% Inputs:
%   K
%   NewDist
%   Closed
% Outputs:
%   Peaks
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.


[PIndPos,PIndNeg] = FindPeaks(K, NewDist, Closed);
if size(PIndPos,1) >= 5 && size(PIndNeg,1) >=2
    Peaks = [PIndPos(1:5,:);PIndNeg(1:2,:)];
    
    clear PIndPos PIndNeg
    
    [~, Ind] = sort(Peaks(:,2));
    Peaks = Peaks(Ind,:);
    clear Ind
    
    m=1;
    while Peaks(1,3) > 0 || Peaks(3,3) > 0
    Peaks = circshift(Peaks,1,1);
    m=m+1;
    if m==8
        %errordlg('Curvature negative peaks misfit. Breaking out!');
        Peaks = [];
        break
    end
        
    end
else
    Peaks = [];
end
end

