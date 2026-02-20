% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


% RamLineWF.m
% Purpose:
%   Descriptor computation stage. Uses the curves/surfaces extracted by
%   Rak_Line_Extraction_WF.m to compute curvature-, torsion-, and plane-based
%   quantitative variables, and appends them to the per-specimen results
%   structure (Res).
%
% Inputs (expected in the workspace from Rak_Line_Extraction_WF.m):
%   Res, Rotv, RotvRam, RotfRam (and associated intermediate variables).
%
% Outputs:
%   Updates Res with the derived variables and curve-level summaries.
%
% Notes:
%   - This script assumes Rak_Line_Extraction_WF.m has been run successfully in
%     the same MATLAB session.
%   - Computational logic and parameter values are intentionally unchanged.

if size(Res,2) == 1
    Res(2) = Res(1);
    Res(2).Side = "N/A";
    RotvRam(2,1) = RotvRam(1,1);
    RotfRam(2,1) = RotfRam(1,1);
end


%Get line on each side
NtchLnR = Rotv(RotvRam{1}(Res(1).NtchLn),:);
NtchLnL = Rotv(RotvRam{2}(Res(2).NtchLn),:);
NtchLnR = [smooth(NtchLnR(:,1),10,'rlowess'),smooth(NtchLnR(:,2),10,'rlowess'),smooth(NtchLnR(:,3),10,'rlowess')];
NtchLnL = [smooth(NtchLnL(:,1),10,'rlowess'),smooth(NtchLnL(:,2),10,'rlowess'),smooth(NtchLnL(:,3),10,'rlowess')];
NRTotLen = sum(vecnorm(NtchLnR(2:end,:)-NtchLnR(1:end-1,:),2,2));
NLTotLen = sum(vecnorm(NtchLnL(2:end,:)-NtchLnL(1:end-1,:),2,2));
NRDists = cumsum(vecnorm(NtchLnR(2:end,:)-NtchLnR(1:end-1,:),2,2));
NLDists = cumsum(vecnorm(NtchLnL(2:end,:)-NtchLnL(1:end-1,:),2,2));

% This section consists of the analytical approach by which the curve is
% expressed in terms of a Fourier Series coefficients. This approach has
% been replaced with the numeric approach nuerically approximating the
% differentiation of the coordinates.
% % %Fourier seires fitting
% % [NRTotLen, NRCoefsX, NRCoefsY, NRCoefsZ, NRGCoefs] = FourierCoefs(NtchLnR, true);
% % [NLTotLen, NLCoefsX, NLCoefsY, NLCoefsZ, NLGCoefs] = FourierCoefs(NtchLnL, true);
% % 
% % %Curve length parameter
% % NRDists = (0:NRTotLen/500:NRTotLen)';
% % NLDists = (0:NLTotLen/500:NLTotLen)';
% % 
% % %Smoothing factor
% % NRWn = 1./(1+exp(((1:1:size(NRCoefsX))-fix(0.25*size(NRCoefsX,1)))/1))';
% % NLWn = 1./(1+exp(((1:1:size(NLCoefsX))-fix(0.25*size(NLCoefsX,1)))/1))';
% % % NRWn = 1./(1+exp(((1:1:size(NRCoefsX))-25)/5))';
% % % NLWn = 1./(1+exp(((1:1:size(NLCoefsX))-25)/5))';
% % 
% % %Fourier series prediction
% % for i = 1:size(NRDists,1)
% %     NRXpred(i,1) = Sum4Fourier(NRTotLen,NRCoefsX,NRDists(i),NRWn);
% %     NRYpred(i,1) = Sum4Fourier(NRTotLen,NRCoefsY,NRDists(i),NRWn);
% %     NRZpred(i,1) = Sum4Fourier(NRTotLen,NRCoefsZ,NRDists(i),NRWn);
% % 
% %     NLXpred(i,1) = Sum4Fourier(NLTotLen,NLCoefsX,NLDists(i),NLWn);
% %     NLYpred(i,1) = Sum4Fourier(NLTotLen,NLCoefsY,NLDists(i),NLWn);
% %     NLZpred(i,1) = Sum4Fourier(NLTotLen,NLCoefsZ,NLDists(i),NLWn);
% % end
% % 
% % %Polynomial prediction for open curves
% % [NRGXs, NRGYs, NRGZs] = GPred(NRDists, NRGCoefs, 0);
% % [NLGXs, NLGYs, NLGZs] = GPred(NLDists, NLGCoefs, 0);
% % 
% % %Fourier transformed line coordinates
% % NRCoords = [NRXpred+NRGXs, NRYpred+NRGYs, NRZpred+NRGZs];
% % NLCoords = [NLXpred+NLGXs, NLYpred+NLGYs, NLZpred+NLGZs];
% % 
% % %Frenet-Serret variables
% % [NRTan, NRNor, NRBin, NRcurv, NRtor] = FrSr(NRDists,NRWn,NRCoefsX, NRCoefsY, NRCoefsZ, NRTotLen, NRGCoefs, true);
% % [NLTan, NLNor, NLBin, NLcurv, NLtor] = FrSr(NLDists,NLWn,NLCoefsX, NLCoefsY, NLCoefsZ, NLTotLen, NLGCoefs, true);

% [NRTan, NRNor, NRBin, NRcurv, NRtor] = FrSrNum3(NtchLnR);
% [NLTan, NLNor, NLBin, NLcurv, NLtor] = FrSrNum3(NtchLnL);

% Resample notch curves to (approximately) equidistant points.
% NOTE: For curve-level descriptors, we parameterize by arc length computed
% from the final curve coordinates used in the differential-geometry step.
[NtchLnR, ~] = EquiDistPT(NtchLnR,true);
[NtchLnL, ~] = EquiDistPT(NtchLnL,true);

% Arc-length parameterization (from the processed curves)
NRDists = [0; cumsum(vecnorm(diff(NtchLnR),2,2))];
NLDists = [0; cumsum(vecnorm(diff(NtchLnL),2,2))];
NRTotLen = NRDists(end);
NLTotLen = NLDists(end);

hNR = mean(vecnorm(diff(NtchLnR),2,2));
hNL = mean(vecnorm(diff(NtchLnL),2,2));
[NRTan, NRNor, NRBin, NRcurv, NRtor] = TNB(NtchLnR(:,1),NtchLnR(:,2),NtchLnR(:,3),hNR);
[NLTan, NLNor, NLBin, NLcurv, NLtor] = TNB(NtchLnL(:,1),NtchLnL(:,2),NtchLnL(:,3),hNL);

% %Smooth curvature and torsion
% NRcurv = filloutliers(NRcurv,"linear","movmedian",7);
% NRcurv = smooth(NRcurv,7,'rlowess');
% NLcurv = filloutliers(NLcurv,"linear","movmedian",7);
% NLcurv = smooth(NLcurv,7,'rlowess');
% NRtor = filloutliers(NRtor,"linear","movmedian",7);
% NRtor = smooth(NRtor,7,'rlowess');
% NLtor = filloutliers(NLtor,"linear","movmedian",7);
% NLtor = smooth(NLtor,7,'rlowess');


NRTanM = mean(NRTan);
NLTanM = mean(NLTan);
NRRamSymDev = abs(acos(dot(NRTanM,Res(1).RamSymVec',2))-0.5*pi);
NLRamSymDev = abs(acos(dot(NLTanM,Res(2).RamSymVec',2))-0.5*pi);


%Integrate curvature and torsion
% NRcurvInt = sum(NRcurv(2:end-1)./(2.*NRDists(2)));
% NLcurvInt = sum(NLcurv(2:end-1)./(2.*NLDists(2)));

% Length-normalized integrals of curvature and torsion along arc length.
% Definitions (for a curve C(s), s in [0,L]):
%   MeanCurvature   = (1/L) * ∫ k(s) ds
%   MeanTorsion     = (1/L) * ∫ τ(s) ds
%   MeanAbsTorsion  = (1/L) * ∫ |τ(s)| ds
%
% Numerical implementation uses the arc-length parameterization computed
% from the processed curve coordinates used in TNB.

NRcurvInt    = trapz(NRDists, NRcurv) ./ NRTotLen;
NLcurvInt    = trapz(NLDists, NLcurv) ./ NLTotLen;

NRtorInt     = trapz(NRDists, NRtor) ./ NRTotLen;
NLtorInt     = trapz(NLDists, NLtor) ./ NLTotLen;

NRAbstorInt  = trapz(NRDists, abs(NRtor)) ./ NRTotLen;
NLAbstorInt  = trapz(NLDists, abs(NLtor)) ./ NLTotLen;


%find maximum central curvature point

NRcurvPeaks = [NRDists(islocalmax(NRcurv))./NRTotLen, NRcurv(islocalmax(NRcurv))];
NRcurvPeaks = NRcurvPeaks(NRcurvPeaks(:,1)>0.2 & NRcurvPeaks(:,1)<0.9,:);
NRcurvPLoc = NRcurvPeaks(NRcurvPeaks(:,2)==max(NRcurvPeaks(:,2)),1);
NRcurvPVal = NRcurvPeaks(NRcurvPeaks(:,2)==max(NRcurvPeaks(:,2)),2);

NLcurvPeaks = [NLDists(islocalmax(NLcurv))./NLTotLen, NLcurv(islocalmax(NLcurv))];
NLcurvPeaks = NLcurvPeaks(NLcurvPeaks(:,1)>0.2 & NLcurvPeaks(:,1)<0.9,:);
NLcurvPLoc = NLcurvPeaks(NLcurvPeaks(:,2)==max(NLcurvPeaks(:,2)),1);
NLcurvPVal = NLcurvPeaks(NLcurvPeaks(:,2)==max(NLcurvPeaks(:,2)),2);

%Cond Line procedure 

CondLnR = Rotv(RotvRam{1}(Res(1).CondLn),:);
CondLnL = Rotv(RotvRam{2}(Res(2).CondLn),:);
CondLnR = [smooth(CondLnR(:,1),7,'rlowess'),smooth(CondLnR(:,2),7,'rlowess'),smooth(CondLnR(:,3),7,'rlowess')];
CondLnL = [smooth(CondLnL(:,1),7,'rlowess'),smooth(CondLnL(:,2),7,'rlowess'),smooth(CondLnL(:,3),7,'rlowess')];
CRTotLen = sum(vecnorm(CondLnR(2:end,:)-CondLnR(1:end-1,:),2,2));
CLTotLen = sum(vecnorm(CondLnL(2:end,:)-CondLnL(1:end-1,:),2,2));
CRDists = cumsum(vecnorm(CondLnR(2:end,:)-CondLnR(1:end-1,:),2,2));
CLDists = cumsum(vecnorm(CondLnL(2:end,:)-CondLnL(1:end-1,:),2,2));

% %Fourier seires fitting
% [CRTotLen, CRCoefsX, CRCoefsY, CRCoefsZ, ~] = FourierCoefs(CondLnR, false);
% [CLTotLen, CLCoefsX, CLCoefsY, CLCoefsZ, ~] = FourierCoefs(CondLnL, false);
% 
% 
% %Curve length parameter
% CRDists = (0:CRTotLen/500:CRTotLen)';
% CLDists = (0:CLTotLen/500:CLTotLen)';
% 
% %Smoothing factor
% CRWn = 1./(1+exp(((1:1:size(CRCoefsX))-fix(0.25*size(CRCoefsX,1)))/1))';
% CLWn = 1./(1+exp(((1:1:size(CLCoefsX))-fix(0.25*size(CLCoefsX,1)))/1))';
% % CRWn = 1./(1+exp(((1:1:size(CRCoefsX))-25)/5))';
% % CLWn = 1./(1+exp(((1:1:size(CLCoefsX))-25)/5))';
% 
% %Fourier series prediction
% for i = 1:size(CRDists,1)
%     CRXpred(i,1) = Sum4Fourier(CRTotLen,CRCoefsX,CRDists(i),CRWn);
%     CRYpred(i,1) = Sum4Fourier(CRTotLen,CRCoefsY,CRDists(i),CRWn);
%     CRZpred(i,1) = Sum4Fourier(CRTotLen,CRCoefsZ,CRDists(i),CRWn);
% 
%     CLXpred(i,1) = Sum4Fourier(CLTotLen,CLCoefsX,CLDists(i),CLWn);
%     CLYpred(i,1) = Sum4Fourier(CLTotLen,CLCoefsY,CLDists(i),CLWn);
%     CLZpred(i,1) = Sum4Fourier(CLTotLen,CLCoefsZ,CLDists(i),CLWn);
% end
% CRCoords = [CRXpred, CRYpred, CRZpred];
% CLCoords = [CLXpred, CLYpred, CLZpred];
% 
% 
% [CRTan, CRNor, CRBin, CRcurv, CRtor] = FrSr(CRDists,CRWn,CRCoefsX, CRCoefsY, CRCoefsZ, CRTotLen, [], false);
% [CLTan, CLNor, CLBin, CLcurv, CLtor] = FrSr(CLDists,CLWn,CLCoefsX, CLCoefsY, CLCoefsZ, CLTotLen, [], false);
% [CRTan, CRNor, CRBin, CRcurv, CRtor] = FrSrNum3(CondLnR);
% [CLTan, CLNor, CLBin, CLcurv, CLtor] = FrSrNum3(CondLnL);
[CondLnR, ~] = EquiDistPT(CondLnR,false);
[CondLnL, ~] = EquiDistPT(CondLnL,false);

% Arc-length parameterization (from the processed curves)
CRDists = [0; cumsum(vecnorm(diff(CondLnR),2,2))];
CLDists = [0; cumsum(vecnorm(diff(CondLnL),2,2))];
CRTotLen = CRDists(end);
CLTotLen = CLDists(end);

hCR = mean(vecnorm(diff(CondLnR),2,2));
hCL = mean(vecnorm(diff(CondLnL),2,2));
[CRTan, CRNor, CRBin, CRcurv, CRtor] = TNB(CondLnR(:,1),CondLnR(:,2),CondLnR(:,3),hCR);
[CLTan, CLNor, CLBin, CLcurv, CLtor] = TNB(CondLnL(:,1),CondLnL(:,2),CondLnL(:,3),hCL);

%Smooth curvature and torsion
CRcurv = filloutliers(CRcurv,"linear","movmedian",7);
CRcurv = smooth(CRcurv,7,'rlowess');
CLcurv = filloutliers(CLcurv,"linear","movmedian",7);
CLcurv = smooth(CLcurv,7,'rlowess');
CRtor = filloutliers(CRtor,"linear","movmedian",7);
CRtor = smooth(CRtor,7,'rlowess');
CLtor = filloutliers(CLtor,"linear","movmedian",7);
CLtor = smooth(CLtor,7,'rlowess');


CRcurvInt    = trapz(CRDists, CRcurv) ./ CRTotLen;
CLcurvInt    = trapz(CLDists, CLcurv) ./ CLTotLen;

CRtorInt     = trapz(CRDists, CRtor) ./ CRTotLen;
CLtorInt     = trapz(CLDists, CLtor) ./ CLTotLen;

CRAbstorInt  = trapz(CRDists, abs(CRtor)) ./ CRTotLen;
CLAbstorInt  = trapz(CLDists, abs(CLtor)) ./ CLTotLen;



OCRCoords = CondLnR - mean(CondLnR);
CRCoef = [OCRCoords(:,1:2), ones(size(OCRCoords,1),1)]\OCRCoords(:,3); % z = a*x + b*y + c
CRPrOZ = OCRCoords(:,1).*CRCoef(1) + OCRCoords(:,2).*CRCoef(2) + CRCoef(3);
CRPlnRMSE = sqrt(mean((OCRCoords(:,3) - CRPrOZ).^2));
CRn = [-CRCoef(1); -CRCoef(2); 1]; % plane normal for z = a*x + b*y + c
CRn = CRn./norm(CRn);
CondSupVecR = Res(1).CondSupVec(:);
CondSupVecR = CondSupVecR./norm(CondSupVecR);
CRPlnAng = acos(max(min(dot(CRn,CondSupVecR),1),-1));

OCLCoords = CondLnL - mean(CondLnL);
CLCoef = [OCLCoords(:,1:2), ones(size(OCLCoords,1),1)]\OCLCoords(:,3); % z = a*x + b*y + c
CLPrOZ = OCLCoords(:,1).*CLCoef(1) + OCLCoords(:,2).*CLCoef(2) + CLCoef(3);
CLPlnRMSE = sqrt(mean((OCLCoords(:,3) - CLPrOZ).^2));
CLn = [-CLCoef(1); -CLCoef(2); 1]; % plane normal for z = a*x + b*y + c
CLn = CLn./norm(CLn);
CondSupVecL = Res(2).CondSupVec(:);
CondSupVecL = CondSupVecL./norm(CondSupVecL);
CLPlnAng = acos(max(min(dot(CLn,CondSupVecL),1),-1));



Res(1).NotchRamSymDev = NRRamSymDev;
Res(1).NotchCurvInt = NRcurvInt;
Res(1).NotchTorInt = NRtorInt;
Res(1).NotchAbsTorInt = NRAbstorInt;
Res(1).NotchPeakCurvLoc = NRcurvPLoc;
Res(1).NotchPeakCurvVal = NRcurvPVal;
Res(1).CondCurvInt = CRcurvInt;
Res(1).CondTorInt = CRtorInt;
Res(1).CondAbsTorInt = CRAbstorInt;
Res(1).CondPlnRMSE = CRPlnRMSE;
Res(1).CondPlnAng = CRPlnAng;

Res(2).NotchRamSymDev = NLRamSymDev;
Res(2).NotchCurvInt = NLcurvInt;
Res(2).NotchTorInt = NLtorInt;
Res(2).NotchAbsTorInt = NLAbstorInt;
Res(2).NotchPeakCurvLoc = NLcurvPLoc;
Res(2).NotchPeakCurvVal = NLcurvPVal;
Res(2).CondCurvInt = CLcurvInt;
Res(2).CondTorInt = CLtorInt;
Res(2).CondAbsTorInt = CLAbstorInt;
Res(2).CondPlnRMSE = CLPlnRMSE;
Res(2).CondPlnAng = CLPlnAng;

if ~exist('MandiblesWF_NoFigures','var') || ~MandiblesWF_NoFigures

Fig1 = figure();
Ax1 = axes(Fig1);
hold(Ax1, 'on');
MadnPlt = patch(Ax1, 'Faces', f, 'Vertices', Rotv,...
    'FaceColor', [1 1 1], 'LineStyle','none','AmbientStrength',0.1,...
    'SpecularExponent',30,'SpecularStrength',0.1, 'FaceLighting', 'gouraud');
L1 = light(Ax1,'Position',[0,40,-60]);
L2 = light(Ax1,'Position',[0,-40,-60]);
L3 = light(Ax1,'Position',[0,40,60]);
L4 = light(Ax1,'Position',[0,-40,60]);
axis(Ax1, 'equal');

Ram3R1 = plot3(Ax1,NtchLnR(:,1),NtchLnR(:,2),NtchLnR(:,3),'LineWidth',1.5,'Color','r');
Ram3L1 = plot3(Ax1,NtchLnL(:,1),NtchLnL(:,2),NtchLnL(:,3),'LineWidth',1.5,'Color','r');
Ram3RC1 = scatter3(Ax1,NtchLnR(:,1),NtchLnR(:,2),NtchLnR(:,3),25,NRcurv,'o','filled');
Ram3LC1 = scatter3(Ax1,NtchLnL(:,1),NtchLnL(:,2),NtchLnL(:,3),25,NLcurv,'o','filled');
Cond3R1 = plot3(Ax1, CondLnR(:,1),CondLnR(:,2),CondLnR(:,3), 'LineWidth',1.5,'Color','r');
Cond3L1 = plot3(Ax1, CondLnL(:,1),CondLnL(:,2),CondLnL(:,3), 'LineWidth',1.5,'Color','r');
Cond3RC1 = scatter3(Ax1,CondLnR(:,1),CondLnR(:,2),CondLnR(:,3),25,CRcurv,'o','filled');
Cond3LC1 = scatter3(Ax1,CondLnL(:,1),CondLnL(:,2),CondLnL(:,3),25,CLcurv,'o','filled');
Tit1 = title(Ax1,'Curvature');
Col1 = colorbar(Ax1);

Fig2 = figure();
Ax2 = axes(Fig2);
hold(Ax2, 'on');
MadnPlt = patch(Ax2, 'Faces', f, 'Vertices', Rotv,...
    'FaceColor', [1 1 1], 'LineStyle','none','AmbientStrength',0.1,...
    'SpecularExponent',30,'SpecularStrength',0.1, 'FaceLighting', 'gouraud');
L1 = light(Ax2,'Position',[0,40,-60]);
L2 = light(Ax2,'Position',[0,-40,-60]);
L3 = light(Ax2,'Position',[0,40,60]);
L4 = light(Ax2,'Position',[0,-40,60]);
axis(Ax2, 'equal');

Ram3R2 = plot3(Ax2,NtchLnR(:,1),NtchLnR(:,2),NtchLnR(:,3),'LineWidth',1.5,'Color','r');
Ram3L2 = plot3(Ax2,NtchLnL(:,1),NtchLnL(:,2),NtchLnL(:,3),'LineWidth',1.5,'Color','r');
Ram3RC2 = scatter3(Ax2,NtchLnR(:,1),NtchLnR(:,2),NtchLnR(:,3),25,NRtor,'o','filled');
Ram3LC2 = scatter3(Ax2,NtchLnL(:,1),NtchLnL(:,2),NtchLnL(:,3),25,NLtor,'o','filled');
Cond3R2 = plot3(Ax2, CondLnR(:,1),CondLnR(:,2),CondLnR(:,3), 'LineWidth',1.5,'Color','r');
Cond3L2 = plot3(Ax2, CondLnL(:,1),CondLnL(:,2),CondLnL(:,3), 'LineWidth',1.5,'Color','r');
Cond3RC2 = scatter3(Ax2,CondLnR(:,1),CondLnR(:,2),CondLnR(:,3),25,CRtor,'o','filled');
Cond3LC2 = scatter3(Ax2,CondLnL(:,1),CondLnL(:,2),CondLnL(:,3),25,CLtor,'o','filled');
Tit2 = title(Ax2,'Torsion');
Col2 = colorbar(Ax2);

end

% -------------------------------------------------------------------------
% Post-run exports (requested):
%   1) Prompt to export extracted curves as red point-only PLY files.
%   2) Prompt to export the full MATLAB workspace as a MAT file.
%
% Notes:
%   - Curve exports are written in the ORIGINAL coordinate system of the input
%     mandible mesh (v0), so they overlay the raw model in MeshLab.
%   - The prompts are shown after every specimen run, whether invoked
%     interactively or via the batch wrapper.
% -------------------------------------------------------------------------
try
    % Resolve output folder + specimen base name
    if exist('File','var') && exist('Path','var') && ~isempty(File) && ~isempty(Path)
        [~, WF_SpecimenBase, ~] = fileparts(File);
        WF_OutDir = Path;
    elseif exist('MandiblesWF_CurrentFile','var') && ~isempty(MandiblesWF_CurrentFile)
        [WF_OutDir, WF_SpecimenBase, ~] = fileparts(MandiblesWF_CurrentFile);
        WF_OutDir = [WF_OutDir, filesep];
    else
        WF_SpecimenBase = 'specimen';
        WF_OutDir = pwd;
    end

    % Prompt 1: export curves
    WF_ChoiceCurves = questdlg(...
        'Export extracted curves as point-only PLY files (red points)?', ...
        'Export curves', 'Yes', 'No', 'No');

    if strcmp(WF_ChoiceCurves,'Yes')
        % Choose reference coordinates for export:
        % Prefer v0 (original coordinates) so the curve overlays the raw mesh.
        if exist('v0','var') && ~isempty(v0)
            WF_Vref = v0;
        else
            WF_Vref = v; % fallback (centered)
        end
        local_export_curves_point_ply(WF_SpecimenBase, WF_OutDir, Res, RotvRam, WF_Vref);
    end

    % Prompt 2: export full workspace
    WF_ChoiceMAT = questdlg(...
        'Export the full MATLAB workspace to a MAT file?', ...
        'Export workspace', 'Yes', 'No', 'No');

    if strcmp(WF_ChoiceMAT,'Yes')
        WF_MatFile = fullfile(WF_OutDir, [WF_SpecimenBase, '_workspace.mat']);
        save(WF_MatFile, '-v7.3');
    end
catch ME
    warning('RamLineWF:PostRunExportFailed', ...
            'Post-run export step failed: %s', ME.message);
end

% CorPts = scatter3(Ax1,[RCoords(RCorPt(1),1);LCoords(LCorPt(1),1)],[RCoords(RCorPt(1),2);LCoords(LCorPt(1),2)],...
%     [RCoords(RCorPt(1),3);LCoords(LCorPt(1),3)],35,'pb','filled');
% NotPts = scatter3(Ax1,[RCoords(RNotPt(1),1);LCoords(LNotPt(1),1)],[RCoords(RNotPt(1),2);LCoords(LNotPt(1),2)],...
%     [RCoords(RNotPt(1),3);LCoords(LNotPt(1),3)],35,'pb','filled');
% % ConPts = scatter3(Ax1,[RXpred(RConPt(1));LXpred(LConPt(1))],[RYpred(RConPt(1));LYpred(LConPt(1))],...
% %     [RZpred(RConPt(1));LZpred(LConPt(1))],35,'pb','filled');


% % Condyle line projection
% %Nearest triangles
% NNs = knnsearch(Rotv,CRCoords,'K',125);
% for i=1:size(NNs,1)
%     CRNNei{i,1} = unique(cell2mat(vertexAttachments(Tr,NNs(i,:)')'))';
% end
% %Add the Condline points in first column
% CRNNei = [mat2cell(CRCoords,ones(size(CRCoords,1),1),3),CRNNei];
% %Their incenters
% CRNNei(:,3) = cellfun(@(x) incenter(Tr,x), CRNNei(:,2),'UniformOutput',false);
% %Their Normals
% CRNNei(:,4) = cellfun(@(x) faceNormal(Tr,x), CRNNei(:,2),'UniformOutput',false);
% %Condylar line point projection onto the plane of each face
% CRNNei(:,5) = cellfun(@(x,y) x-y,CRNNei(:,1),CRNNei(:,3),'UniformOutput',false);
% CRNNei(:,5) = cellfun(@(x,y) dot(x,y,2),CRNNei(:,5),CRNNei(:,4),'UniformOutput',false);
% CRNNei(:,5) = cellfun(@(x,y,z) x-y.*z,CRNNei(:,1),CRNNei(:,5),CRNNei(:,4),'UniformOutput',false);
% %Check if projection in triangle
% CRNNei(:,6) = cellfun(@(x) [Rotv(f(x,1),:),Rotv(f(x,2),:),Rotv(f(x,3),:)],CRNNei(:,2),'UniformOutput',false);
% 
% %In place X,Y,Z, all
% CRNNei(:,7) = cellfun(@(x,y) all([any(x(:,1)<[y(:,1),y(:,4),y(:,7)],2),any(x(:,1)>[y(:,1),y(:,4),y(:,7)],2)],2),CRNNei(:,5),CRNNei(:,6),'UniformOutput',false);
% CRNNei(:,8) = cellfun(@(x,y) all([any(x(:,2)<[y(:,2),y(:,5),y(:,8)],2),any(x(:,2)>[y(:,2),y(:,5),y(:,8)],2)],2),CRNNei(:,5),CRNNei(:,6),'UniformOutput',false);
% CRNNei(:,9) = cellfun(@(x,y) all([any(x(:,3)<[y(:,3),y(:,6),y(:,9)],2),any(x(:,3)>[y(:,3),y(:,6),y(:,9)],2)],2),CRNNei(:,5),CRNNei(:,6),'UniformOutput',false);
% CRNNei(:,10) = cellfun(@(x,y,z) all([x,y,z],2),CRNNei(:,7),CRNNei(:,8),CRNNei(:,9),'UniformOutput',false); 
% 
% CRNNei = [CRNNei(:,1), cellfun(@(x,y) x(y,:),CRNNei(:,5),CRNNei(:,10),'UniformOutput',false)];
% CRNNei(:,3) = cellfun(@(x,y) pdist2(x,y),CRNNei(:,1),CRNNei(:,2),'UniformOutput',false); 
% CRNNei(:,3) = cellfun(@(x) find(x==min(x),1),CRNNei(:,3),'UniformOutput',false);
% 
% CRCoordsPrj = zeros(size(CRCoords,1),3);
% for i=1:size(CRNNei,1)
%     if ~isempty(CRNNei{i,3})
%         CRCoordsPrj(i,:) = CRNNei{i,2}(CRNNei{i,3},:);
%     else
%         CRCoordsPrj(i,:) = CRNNei{i,1};
%     end
% end


% -------------------------------------------------------------------------
% Local helper: export extracted curves (notch + condylar delineation) as
% point-only PLY files with red vertex colors.
%
% Output convention (per available side):
%   <SpecimenBase>_notch_<R/L>_points.ply
%   <SpecimenBase>_condyle_<R/L>_points.ply
% -------------------------------------------------------------------------
function local_export_curves_point_ply(specimenBase, outDir, Res, RotvRam, Vref)

% Validate required variables from the workflow
if nargin < 5 || isempty(Res) || isempty(RotvRam) || isempty(Vref)
    error('Cannot export curves: missing required inputs (Res / RotvRam / Vref).')
end

% Determine which Res entries are real sides (exclude synthetic N/A row)
if isfield(Res,'Side')
    sides = string({Res.Side});
    valid = sides ~= "N/A";
else
    sides = repmat("", 1, numel(Res));
    valid = true(1, numel(Res));
end

for i = 1:numel(Res)
    if ~valid(i)
        continue
    end

    % Side suffix for filenames
    sideSuffix = "";
    if strlength(sides(i)) > 0
        if sides(i) == "Right"
            sideSuffix = "R";
        elseif sides(i) == "Left"
            sideSuffix = "L";
        else
            % sanitize arbitrary side strings
            sideSuffix = regexprep(char(sides(i)), '\\s+', '');
        end
    end

    % Export notch curve (indices into RotvRam{i})
    if isfield(Res,'NtchLn') && ~isempty(Res(i).NtchLn)
        idxGlobal = RotvRam{i}(Res(i).NtchLn(:));
        idxGlobal = idxGlobal(:);
        idxGlobal = idxGlobal(idxGlobal >= 1 & idxGlobal <= size(Vref,1));
        idxGlobal = unique(idxGlobal, 'stable');
        pts = Vref(idxGlobal,:);

        if ~isempty(sideSuffix)
            outFile = fullfile(outDir, sprintf('%s_notch_%s_points.ply', specimenBase, sideSuffix));
        else
            outFile = fullfile(outDir, sprintf('%s_notch_points.ply', specimenBase));
        end

        local_write_point_ply_red(outFile, pts);
    end

    % Export condylar delineation curve
    if isfield(Res,'CondLn') && ~isempty(Res(i).CondLn)
        idxGlobal = RotvRam{i}(Res(i).CondLn(:));
        idxGlobal = idxGlobal(:);
        idxGlobal = idxGlobal(idxGlobal >= 1 & idxGlobal <= size(Vref,1));
        idxGlobal = unique(idxGlobal, 'stable');
        pts = Vref(idxGlobal,:);

        if ~isempty(sideSuffix)
            outFile = fullfile(outDir, sprintf('%s_condyle_%s_points.ply', specimenBase, sideSuffix));
        else
            outFile = fullfile(outDir, sprintf('%s_condyle_points.ply', specimenBase));
        end

        local_write_point_ply_red(outFile, pts);
    end
end

end


% -------------------------------------------------------------------------
% Local helper: write a point-only PLY with per-vertex red colors.
%   - No faces (face count = 0)
%   - Coordinates written in the provided reference frame
%
% Preference order:
%   (1) surfaceMesh + writeSurfaceMesh (if available)
%   (2) pointCloud + pcwrite (if available)
%   (3) manual ASCII PLY writer (fallback)
% -------------------------------------------------------------------------
function local_write_point_ply_red(outFile, pts)

if isempty(pts)
    return
end

pts = double(pts);
col = uint8(repmat([255 0 0], size(pts,1), 1));

% Attempt 1: surfaceMesh ecosystem
try
    if exist('writeSurfaceMesh','file') == 2
        try
            % Try common constructor signatures
            try
                SMp = surfaceMesh(zeros(0,3), pts); % (Faces, Vertices)
            catch
                SMp = surfaceMesh(pts, zeros(0,3)); % (Vertices, Faces)
            end

            % Attach vertex colors where supported
            if isprop(SMp,'VertexColors')
                SMp.VertexColors = col;
            end

            writeSurfaceMesh(SMp, outFile);
            return
        catch
            % Fall through to the next method
        end
    end
catch
end

% Attempt 2: pointCloud ecosystem
try
    if exist('pointCloud','class') == 8 && exist('pcwrite','file') == 2
        pc = pointCloud(pts, 'Color', col);
        pcwrite(pc, outFile, 'Encoding', 'ascii');
        return
    end
catch
end

% Attempt 3: manual ASCII PLY
fid = fopen(outFile, 'w');
if fid < 0
    error('Cannot open file for writing: %s', outFile)
end

fprintf(fid, 'ply\n');
fprintf(fid, 'format ascii 1.0\n');
fprintf(fid, 'element vertex %d\n', size(pts,1));
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'property uchar red\n');
fprintf(fid, 'property uchar green\n');
fprintf(fid, 'property uchar blue\n');
fprintf(fid, 'element face 0\n');
fprintf(fid, 'property list uchar int vertex_indices\n');
fprintf(fid, 'end_header\n');

for k = 1:size(pts,1)
    fprintf(fid, '%.6f %.6f %.6f %d %d %d\n', ...
        pts(k,1), pts(k,2), pts(k,3), col(k,1), col(k,2), col(k,3));
end

fclose(fid);

end
