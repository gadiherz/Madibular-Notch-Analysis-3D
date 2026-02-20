% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Rotv] = PosOnTable(v,f)
% PosOnTable is a core function of the workflow used to rotate the
% mandibale about the X axis normal to the Symmetry plane so that its base
% would be positioned on a "table" perpendicular to SymPlane
%
% Inputs:
%   v
%   f
% Outputs:
%   Rotv
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%   The function takes as input the FullDB after SymmTens_2 and adds to it
%   the following fields: rotation angle (FinAng), rotated vertices (Rotv),
%   the positioning option selected (PosSel), whether it was flipped about
%   Y axis (flipped). It adds three structs (Neg, Pos, Comb) contating the
%   curve data of each projection (coords, length, peaks) the base angle
%   and points through which CommTan, and ramus angle and points. If Pos
%   and Neg agree (5 deg diff) the process is automatic and Comb is empty.
%   If not, the user selects the best option and Comb is created. In this
%   situation it is possible that either of the option structs will not
%   contatin all fields in case no good configuration was found.

      
%negative side - gets points of negative half

%takes the positive half of X (i.e. half positive side of symmetry
%plane). Runs Fourier transform on its ortho projection onto the
%Symplane, first with high smoothing parameter to find approximate
%peaks, than with lower smoothing to find the percise peaks. If no good
%configuration is found (p,n,p,n,p,p,p) returns all zeros.
v4Pos = v(v(:,3)<=0,:);
Sz = size(v4Pos,1);
if Sz>40e3
    % Deterministic subsampling (replaces randi-based sampling for reproducibility).
    nSmp = 40e3;
    SmpInd = unique(round(linspace(1,Sz,nSmp))','stable');
else
    SmpInd = (1:Sz)';
end
b4Pos = boundary(v4Pos(SmpInd,2),v4Pos(SmpInd,1));
b4Pos = [v4Pos(SmpInd(b4Pos),2),v4Pos(SmpInd(b4Pos),1)];
[NewDist,~,~,~,~,K,~,~,~] = FourierCoords(b4Pos,6);
PeaksNeg = JustPeaks(K,NewDist,true);
[NewDist,~,~,~,Theta,K,Xpred,Ypred,LenNeg] = FourierCoords(b4Pos,45);
if isempty(PeaksNeg)
    angNeg = 0;
    RamangNeg = 0;
    RotPCNeg = [Xpred,Ypred];
    BotPtsNeg = [0,0;0,0];
    RamPtsNeg = [0,0;0,0];
% finds common tangent of bottom section (always between peaks 5&6 but
% may be from different directions) this gives rotation ang of mandible
% base to table, and the two points (or rearly one) between which the
% com tan passes
else
    
    PeaksNeg = PercPeaks(PeaksNeg,K,NewDist); 
    [angNeg,BotPtsNeg] = FindComTanF2(PeaksNeg(5:6,:),NewDist, Theta, K,Xpred, Ypred);  
    RotPCNeg = [Xpred,Ypred];
    % Checks if angle rotates base to "ceiling" than flip 180 to table
    RotPCNegTmp = ([cos(-angNeg), -sin(-angNeg); sin(-angNeg), cos(-angNeg)] *[Xpred,Ypred]')';
    if RotPCNegTmp(PeaksNeg(5,1),2) > 0 && RotPCNegTmp(PeaksNeg(6,1),2) > 0
        angNeg = angNeg + pi;
        RotPCNegTmp = ([cos(-pi), -sin(-pi); sin(-pi), cos(-pi)] *RotPCNegTmp')';
    end
    %Checks direction: if the first peak has higher Y value than third
    %both negative than the chin is towards negative X and ramus
    %section is between the 6 and 7 peak, otherwise the ramus is
    %between the 4 and 5 peak and chin is towards positive X
    if RotPCNegTmp(PeaksNeg(1,1),2) > RotPCNegTmp(PeaksNeg(3,1),2)
        [RamangNeg,RamPtsNeg] = FindComTanF2(PeaksNeg(6:7,:),NewDist,Theta,K,Xpred,Ypred);
    else
        [RamangNeg,RamPtsNeg] = FindComTanF2(PeaksNeg(4:5,:),NewDist,Theta,K,Xpred,Ypred);

    end
end
clear v4Pos b4Pos Sz SmpInd NewDist Wn CoefsX CoefsY Xs Ys Theta K Xpred Ypred RotPCNegTmp

%positive side - exactly the same process as above - see comments
v4Pos = v(v(:,3)>=0,:);
Sz = size(v4Pos,1);
if Sz>40e3
    % Deterministic subsampling (replaces randi-based sampling for reproducibility).
    nSmp = 40e3;
    SmpInd = unique(round(linspace(1,Sz,nSmp))','stable');
else
    SmpInd = (1:Sz)';
end
b4Pos = boundary(v4Pos(SmpInd,2),v4Pos(SmpInd,1));
b4Pos = [v4Pos(SmpInd(b4Pos),2),v4Pos(SmpInd(b4Pos),1)];
[NewDist,~,~,~,~,K,~,~,~] = FourierCoords(b4Pos,6);
PeaksPos = JustPeaks(K,NewDist,true);
[NewDist,~,~,~,Theta,K,Xpred,Ypred,LenPos] = FourierCoords(b4Pos,45);
if isempty(PeaksPos)
    angPos = 0;
    RamangPos = 0;
    RotPCPos = [Xpred,Ypred];
    BotPtsPos = [0,0;0,0];
    RamPtsPos = [0,0;0,0];
else
    
    PeaksPos = PercPeaks(PeaksPos,K,NewDist);
    [angPos, BotPtsPos]  = FindComTanF2(PeaksPos(5:6,:),NewDist,Theta,K,Xpred,Ypred);
    RotPCPos = [Xpred,Ypred];
    RotPCPosTmp = ([cos(-angPos), -sin(-angPos); sin(-angPos), cos(-angPos)] *[Xpred,Ypred]')';
    if RotPCPosTmp(PeaksPos(5,1),2) > 0 && RotPCPosTmp(PeaksPos(6,1),2) > 0
        angPos = angPos + pi;
        RotPCPosTmp = ([cos(-pi), -sin(-pi); sin(-pi), cos(-pi)] *RotPCPosTmp')';
    end
    if RotPCPosTmp(PeaksPos(1,1),2) > RotPCPosTmp(PeaksPos(3,1),2)
        [RamangPos,RamPtsPos] = FindComTanF2(PeaksPos(6:7,:),NewDist,Theta,K,Xpred,Ypred);
    else
        [RamangPos,RamPtsPos] = FindComTanF2(PeaksPos(4:5,:),NewDist,Theta, K,Xpred,Ypred);
    end
end
clear v4Pos b4Pos Sz SmpInd NewDist Wn CoefsX CoefsY Xs Ys Theta K Xpred Ypred RotPCPosTmp

% for the actual positioning, checks if the base angles of the two
% sides agree within 5 degrees if so the reotation angle for the entire
% mandible is their mean (remember this difference for display later)
if rad2deg(abs(angNeg-angPos))<=5 && angNeg ~= 0 && angPos ~= 0 
    FinAng = mean([angNeg;angPos]);
    Rotv = ([cos(FinAng+0.5*pi), -sin(FinAng+0.5*pi),0;...
        sin(FinAng+0.5*pi), cos(FinAng+0.5*pi),0; 0,0,1]* v')';
% Checks if in 3D chin is towards negative Z if so, it flips the
% mandible 180 degs about the Y axis so chin is towards positive and
% keeps a note in FullDB that it was flipped. This is important because
% if flipped: Neg => Left & Pos => Right
% if not fliped: Neg => Right & Pos => Left
    TriCrds = zeros(size(f,1),3,2);
    TriCrds(:,:,1) = reshape(Rotv(f(:),3),[],3);
    TriCrds(:,:,2) = reshape(Rotv(f(:),1),[],3);
    Tri =  find(all([any(TriCrds(:,:,1)<0,2),any(TriCrds(:,:,1)>0,2)],2),1);
    if all(TriCrds(Tri,:,2)>0)
        Rotv = ([cos(pi), 0, sin(pi); 0, 1, 0; -sin(pi), 0, cos(pi)] *Rotv')';
    end
    clear Tri TriCrds

%If Pos and Neg angles disagree at 5 degrees, runs the same procedure
%on the ortho projection of all vertices onto the SymPlane. Exactly the
%same procedure as Pos and Neg, see comments above.
else
    v4Pos = v;
    Sz = size(v4Pos,1);
    if Sz>40e3
        % Deterministic subsampling (replaces randi-based sampling for reproducibility).
        nSmp = 40e3;
        SmpInd = unique(round(linspace(1,Sz,nSmp))','stable');
    else
        SmpInd = (1:Sz)';
    end
    b4Pos = boundary(v4Pos(SmpInd,2),v4Pos(SmpInd,1));
    b4Pos = [v4Pos(SmpInd(b4Pos),2),v4Pos(SmpInd(b4Pos),1)];
    [NewDist,~,~,~,~,K,~,~,~] = FourierCoords(b4Pos,6);
    PeaksComb = JustPeaks(K,NewDist,true);
    [NewDist,~,~,~,Theta,K,Xpred,Ypred,LenComb] = FourierCoords(b4Pos,45);
   
    if isempty(PeaksComb)
        angComb = 0;
        RamangComb = 0;
        RotPCComb = [Xpred,Ypred];
        BotPtsComb = [0,0;0,0];
        RamPtsComb = [0,0;0,0];
    else
         PeaksComb = PercPeaks(PeaksComb,K,NewDist);
        [angComb, BotPtsComb] = FindComTanF2(PeaksComb(5:6,:),NewDist,Theta,K,Xpred,Ypred);
        RotPCComb = [Xpred,Ypred];
        RotPCCombTmp = ([cos(-angComb), -sin(-angComb); sin(-angComb), cos(-angComb)] *[Xpred,Ypred]')';
        if RotPCCombTmp(PeaksComb(5,1),2) > 0 && RotPCCombTmp(PeaksComb(6,1),2) > 0
            angComb = angComb + pi; 
            RotPCCombTmp = ([cos(-pi), -sin(-pi); sin(-pi), cos(-pi)] *RotPCCombTmp')';
        end
        if RotPCCombTmp(PeaksComb(1,1),2) > RotPCCombTmp(PeaksComb(3,1),2)
            [RamangComb,RamPtsComb] = FindComTanF2(PeaksComb(6:7,:),NewDist,Theta,K,Xpred,Ypred);
        else
            [RamangComb,RamPtsComb] = FindComTanF2(PeaksComb(4:5,:),NewDist,Theta,K,Xpred,Ypred);
        end
    end
    clear NewDist Wn CoefsX CoefsY Xs Ys Theta K Xpred Ypred RotPCCombTmp
    
    %Now applies all three options for user selection. If all options
    %dont provide any good configuration model is skipped.
    if all([angNeg, angPos, angComb]==0)
        h = errordlg({'No good configuration of curvature peaks found.';'Cannot position model on table.'});
        waitfor(h);
        delete(h);
        clear h
    % If there is at least one good configuration user chooses from
    % Select Side GUI. As a rule, if it is possible it's better top
    % choose either Pos or Neg. The actual rotation angle about X is
    % determined according to selection. PosSel indicates which side
    % was selected. Canceling returns 0.
    else
        SideSel = SelectSide4Pos(RotPCNeg, PeaksNeg, angNeg, RotPCPos, PeaksPos, angPos, RotPCComb, PeaksComb, angComb);
        switch SideSel
            case 1
               FinAng = angNeg;
            case 2
                FinAng = angComb;
            case 3
                FinAng = angPos;
            case 4
                FinAng = 0;
        end
        Rotv = ([cos(FinAng+0.5*pi), -sin(FinAng+0.5*pi),0;...
            sin(FinAng+0.5*pi), cos(FinAng+0.5*pi),0;0,0,1]* v')';
        % Rotv field in FullDB contains the fully positioned and 
        % rotated vertices accordig to the selected tensor and angle
        
        % Checks directionof chin and flips about Y as necessary,
        % exactly the same as above, see comments.
        TriCrds = zeros(size(f,1),3,2);
        TriCrds(:,:,1) = reshape(Rotv(f(:),3),[],3);
        TriCrds(:,:,2) = reshape(Rotv(f(:),1),[],3);
        Tri =  find(all([any(TriCrds(:,:,1)<0,2),any(TriCrds(:,:,1)>0,2)],2),1);
        if all(TriCrds(Tri,:,2)>0)
            Rotv = ([cos(pi), 0, sin(pi); 0, 1, 0; -sin(pi), 0, cos(pi)] *Rotv')';
        end
        clear Tri TriCrds
        
        
    end
end 

end

