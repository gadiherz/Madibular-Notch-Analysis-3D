% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


% Rak_Line_Extraction_WF.m
% Purpose:
%   Interactive extraction stage. Loads a 3D mandible mesh, orients it to a
%   quasi-anatomical frame, segments usable rami/condyle regions, and extracts
%   anatomically defined curves (mandibular notch ridge and the superior
%   condylar delineation curve) for downstream differential-geometric
%   quantification.
%
% Inputs:
%   None. The mesh file is selected via a file dialog.
%
% Outputs (workspace variables created/updated by the script):
%   Rotv, Rotf, RotvRam, RotfRam, Res (and intermediate variables used by the
%   workflow). See the Supplementary Online Methods (SOM) for definitions.
%
% Interactive steps:
%   - Select an input mesh in the file dialog.
%   - Follow prompts in CompFragSelect (complete mandible vs fragment; side
%     selection; point picking for fragment positioning).
%
% Notes:
%   - This script is designed for interactive use; user choices are part of the
%     workflow and should be recorded when used for publication.
%   - Computational logic and parameter values are intentionally unchanged.

% NOTE: This workflow stage is primarily interactive. When executed under a
% batch wrapper (MandiblesWF_BatchMode = true), do not clear the caller
% workspace (e.g., accumulated results table). Instead, clear only variables
% that do not belong to the batch controller namespace (MandiblesWF_*).
if ~exist('MandiblesWF_BatchMode','var') || ~MandiblesWF_BatchMode
    clear
    clc
else
    vars = who;
    keep = startsWith(vars,'MandiblesWF_');
    clearvars(vars{~keep});
end
%Select file
if exist('MandiblesWF_BatchMode','var') && MandiblesWF_BatchMode && exist('MandiblesWF_CurrentFile','var') && ~isempty(MandiblesWF_CurrentFile)
    [p0, f0, e0] = fileparts(MandiblesWF_CurrentFile);
    File = [f0, e0];
    Path = [p0, filesep];
else
    [File, Path] = uigetfile({'*.ply;*.stl;*.obj','3D mesh files (*.ply, *.stl, *.obj)';'*.*','All files'}, 'Select 3D model', 'MultiSelect', 'off');
end

if isequal(File,0)
    if ~exist('MandiblesWF_BatchMode','var') || ~MandiblesWF_BatchMode
        clear
    end
    return
end

%Read mesh
tic
%[v,f] = read_ply([Path,File]);
SM = readSurfaceMesh([Path,File]);
v0 = SM.Vertices; % Original coordinates (for curve export in original mesh frame)
v = v0;
f = double(SM.Faces);

% Large-mesh warning (user requested): if the input model exceeds 750k faces,
% notify that processing may take a long time.
if size(f,1) > 750000
    try
        h = warndlg(sprintf(['WARNING: The input mesh contains %d faces (>750,000). ',...
            'Processing may take a long time.'], size(f,1)), ...
            'Large mesh warning', 'modal');
        uiwait(h);
    catch
        warning('Large mesh warning: %d faces (>750,000). Processing may take a long time.', size(f,1));
    end
end

v = v - mean(v);
'load'
toc

%Reduce model to 100k triangles 50k vertices for principal curvature
%calculation
% if size(v,1) > 300e3
%     P = patch('Faces', f, 'Vertices', v);
%     P = reducepatch(P, 300e3);
%     v = P.vertices;
%     f = P.faces;
%     close force all
%     clear P
%     'reduce'
%     toc
% %     [v,f] = CleanRedTr(v,f);
% %     toc
% %     Tr = triangulation(f,v);
% %     fN = faceNormal(Tr);
% %     toc
% %     Tr = FixManifold(Tr);
% %     toc
% %     Tr = FixNorDir(Tr);
% %     toc
% end

%Translate to origin and fix colinear or redundant faces
% v = v - mean(v);
% [v,f] = CleanRedTr(v,f);
% 
% toc

%Triangulate, fix manifoldness and arrange all face normals coherently
Tr = triangulation(f,v);
fN = faceNormal(Tr);
InCnt = incenter(Tr);
% Tr = FixManifold(Tr);
% Tr = FixNorDir(Tr);
% 
% % Flip faces and retriangulate if needed 
% f = Tr.ConnectivityList;
% v = Tr.Points;
% fN = faceNormal(Tr);
% InCnt = incenter(Tr);
% Nf = NorFlpMand(v,f,fN,InCnt);
% if ~isempty(Nf)
%     f = Nf;
%     Tr = triangulation(f,v);
%     fN = faceNormal(Tr);
% end
% clear Nf

%find positioning tensor and orthonormal basis for rotation
%Rotate to symmetry palne
fArea = 0.5*vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2);
Tens = Tens2(fN,fArea);
[EgVec, EgVal] = eig(Tens);
if det(EgVec) < 0
    EgVec = -EgVec;
end
Rotv = (EgVec'*v')';
'find symmetry plane'
toc


% %Complete mandible or fargment
% Fig1 = figure();
% Ax1 = axes(Fig1);
% hold(Ax1, 'on');
% MadnPlt = patch('Faces', f, 'Vertices', Rotv,...
% 'FaceColor',[1,1,1], 'LineStyle','none','AmbientStrength',0.1,...
% 'SpecularExponent',30,'SpecularStrength',0.1, 'FaceLighting', 'gouraud');
% L1 = light(Ax1,'Position',[0,40,-60]);
% L2 = light(Ax1,'Position',[0,-40,-60]);
% L3 = light(Ax1,'Position',[0,40,60]);
% L4 = light(Ax1,'Position',[0,-40,60]);
% axis(Ax1, 'equal');

%UI for complete/fragment identification and secondary positioning
App = CompFragSelect(f, Rotv);
while ~App.Done
    pause(0.05)
end
Rotv = App.RotV; 
RotvRam = App.RotvRam;
RotfRam = App.RotfRam;
Side = App.FrgSd;
ChnPt = App.ChinBasePt;
NchMin = App.NchMin;
delete(App)

%retriangulate
Tr = triangulation(f,Rotv);
'user input & triangulate'
toc
%%



%Calculate principal curvature
[Kmax, KDmax, Kmin, KDmin,vN] = SurfPrCurv(Tr);

%Calculate ridgness

% Kmax = LapSmth(Rotv,Kmax,25,5);
% Kmin = LapSmth(Rotv,Kmin,25,5);
Crvns = CalcRdgnss(Rotv,Kmax,Kmin,KDmax);
%Crvns = LapSmth(Rotv,Crvns,20,3);
Crvns = Crvns + abs(min(Crvns));
CrvTh = quantile(Crvns,100);
'calculate curvature'
toc


%% Find coords (of vertices) on mesh of the Rak Line ridge
if iscell(RotvRam)
    NRuns = 2;
else
    NRuns = 1;
    RotvRam = {RotvRam};
    RotfRam = {RotfRam};
    NchMin = {NchMin};
end
for i=1:NRuns
    %Get adjust principal curvature indices to ramus vertices
    RKmax = Kmax(RotvRam{i});
    RKmin = Kmin(RotvRam{i});
    RKDmax = KDmax(RotvRam{i},:);
    RKDmin = KDmin(RotvRam{i},:);
    RCrvns = Crvns(RotvRam{i});
    RCrvTh = quantile(RCrvns,100);
    RvN = vN(RotvRam{i},:);
    RKGauss = RKmax.*RKmin;
    RKGTh = quantile(RKGauss,100);
    
    [NchMinInd,CrndInd] = CrndNchMin(Rotv,RotvRam{i},RCrvns,RKGauss);
    %Triangulate Ramus
    TrR = triangulation(RotfRam{i},Rotv(RotvRam{i},:));
    %Find high curvature points , ridge points and edges between them
     RCrvPts = find(RCrvns>=RCrvTh(60));
    % RCrvPts = find(RKGauss<RKGTh(5) | RKGauss>RKGTh(95));

  % RdgPtsInd = FindRdgPts(TrR,RCrvns, RKDmax, RCrvPts);
    RdgEdgs = RdgConn(TrR,RCrvPts,RKDmin);
    'find edge points & connection'
    toc
    [~,NewRdgEdgs] = EdgeUni(TrR, RdgEdgs, RKDmin);
    RdgPts = CenterRidge(TrR, NewRdgEdgs, RKDmin);
    'center ridge'
    toc
%     CondLn = CondLine(Rotv,RKGauss,RotvRam{i},RotfRam{i},TrR,RvN,RCrvns,CrndInd);
%     CondLn = CondLineSmth(Rotv,RotvRam{i},CondLn);
    [CondLn, CondCntr, CondSupVec, CondSymVec, CondCrdSys] = CondLnPrjC3(Rotv, RotvRam{i}, RotfRam{i}, RCrvPts, RvN, RCrvns, CrndInd, ChnPt,NchMinInd);
    CondLnNN = nearestNeighbor(TrR,CondLn);
    CondLnDst = vecnorm(CondLn-Rotv(RotvRam{i}(CondLnNN),:),2,2);
    CondLnNN = unique(CondLnNN(CondLnDst<0.5),'stable');
    CondLn = [CondLnNN;CondLnNN(1)];
    [CondLn,CondSurf, TriCrvMean, TriCrvMeanSD, TriCrvGaus, TriCrvGausSD] = CondSurfArea(Rotv, RotvRam{i}, RotfRam{i}, TrR, CondLn,CondCntr, RKmax, RKmin);
%     NtchLn = NotchLine(RvN,RKDmin,RCrvns,CrndInd,RdgPts,Rotv,RotvRam{i},CondLn);
    NtchLn = NotchLine3(Rotv, RotvRam{i},CrndInd,NchMinInd,RCrvPts, RCrvns, CondLn, RKDmin, RvN);
    [RamSymVec, RamCrdSys] = RamSymPln(RvN);
    % %Euler angles and Frobenius norm between coordinate systems
    % [RMX,RMY,RMZ,RMFN] = EulerAngMat([1,0,0;0,1,0;0,0,1],RamCrdSys);
    % [CMX,CMY,CMZ,CMFN] = EulerAngMat([1,0,0;0,1,0;0,0,1],CondCrdSys);
    % [CRX,CRY,CRZ,CRFN] = EulerAngMat(CondCrdSys,RamCrdSys);
    RamMandSymAng = acos(dot(RamSymVec',[1,0,0]));
    if RamMandSymAng > 0.5*pi
        RamMandSymAng = acos(dot(RamSymVec',[-1,0,0]));
    end
    CondRamSymAng = acos(dot(RamSymVec',CondSymVec));
    if CondRamSymAng > 0.5*pi
        CondRamSymAng = acos(dot(RamSymVec',-CondSymVec));
    end
%      Lines = [NtchLn(1:end-1)';circshift(CondLn(1:end-1),-find(CondLn==NtchLn(end)))];
     %Lines = ConRdgPts(RdgPts,RKDmax);
     'build lines'
    toc
% 
% %     Plot
% %     fig3=figure();
% %     ax3=axes(fig3);
% %     patch(ax3, 'Faces', RotfRam{i}, 'Vertices', Rotv(RotvRam{i},:),...
% %     'FaceColor', [1,1,1], 'LineStyle','none','AmbientStrength',0.05,...
% %     'SpecularExponent',30,'SpecularStrength',0.05, 'FaceLighting', 'gouraud');
% %     hold on
% %     L1=light(ax3,'Position',[0,40,-60]);
% %     L2=light(ax3,'Position',[0,-40,-60]);
% %     L3=light(ax3,'Position',[0,40,60]);
% %     L4=light(ax3,'Position',[0,-40,60]);
% %     axis(ax3, 'equal');
% %     for j =1:size(Lines,1)
% %     plot3(Rotv(RotvRam{i}(Lines{j}),1),Rotv(RotvRam{i}(Lines{j}),2),Rotv(RotvRam{i}(Lines{j}),3),'LineWidth',2);
% %     end


%Res(i).RdgEdgs = RdgEdgs;
%Res(i).RdgPts = RdgPts;
%Res(i).Lines = Lines;

%vectors and arrays
Res(i).RCrvPts = RCrvPts;
Res(i).NtchLn = NtchLn;
Res(i).CondLn = CondLn;
Res(i).CondCntr = CondCntr;
Res(i).CondSupVec = CondSupVec;
Res(i).CondSymVec = CondSymVec;
Res(i).RamSymVec = RamSymVec;

%Scalars and strings
if ~iscell(Side)
    Res(i).Side = string(Side);
elseif i==1
    Res(i).Side = "Right";
else
    Res(i).Side = "Left";
end
Res(i).CondSurf = CondSurf;
Res(i).TriCrvMean = TriCrvMean;
Res(i).TriCrvMeanSD = TriCrvMeanSD;
Res(i).TriCrvGaus = TriCrvGaus;
Res(i).TriCrvGausSD = TriCrvGausSD;
Res(i).RamMandSymAng = RamMandSymAng;
Res(i).CondRamSymAng = CondRamSymAng;


end

 %Plot

 fig3=figure();
 ax3=axes(fig3);
 patch(ax3, 'Faces', f, 'Vertices', Rotv,...
     'FaceColor', 'interp', 'CData', filloutliers(Crvns,'clip'), 'LineStyle','none','AmbientStrength',0.05,...
     'SpecularExponent',30,'SpecularStrength',0.05, 'FaceLighting', 'gouraud');
 hold on
 L1=light(ax3,'Position',[0,40,-60]);
 L2=light(ax3,'Position',[0,-40,-60]);
 L3=light(ax3,'Position',[0,40,60]);
 L4=light(ax3,'Position',[0,-40,60]);
 axis(ax3, 'equal');
 Tit= title(ax3, File(1:end-4),'Interpreter','none');
% for i=1:NRuns
%     scatter3(Rotv(RotvRam{i}(cell2mat(Res(i).RdgPts(:,1))),1),Rotv(RotvRam{i}(cell2mat(Res(i).RdgPts(:,1))),2),Rotv(RotvRam{i}(cell2mat(Res(i).RdgPts(:,1))),3),'.r');
% for j =1:size(Res(i).Lines,1)
%     plot3(Rotv(RotvRam{i}(Res(i).Lines{j}),1),Rotv(RotvRam{i}(Res(i).Lines{j}),2),Rotv(RotvRam{i}(Res(i).Lines{j}),3),'LineWidth',2);
%  end
% end
scatter3(ax3,Rotv(ChnPt,1),Rotv(ChnPt,2),Rotv(ChnPt,3),100,'or','filled');
for i=1:NRuns
%scatter3(Rotv(RotvRam{i}(cell2mat(Res(i).RdgPts(:,1))),1),Rotv(RotvRam{i}(cell2mat(Res(i).RdgPts(:,1))),2),Rotv(RotvRam{i}(cell2mat(Res(i).RdgPts(:,1))),3),'.r');
    if any(ismember(Res(i).NtchLn,Res(i).CondLn))
        Line = [Res(i).NtchLn(1:end-1)';circshift(Res(i).CondLn(1:end-1),-find(Res(i).CondLn==Res(i).NtchLn(end)))];
        plot3(Rotv(RotvRam{i}(Line),1),Rotv(RotvRam{i}(Line),2),Rotv(RotvRam{i}(Line),3),'r','LineWidth',2);
    else
        plot3(Rotv(RotvRam{i}(Res(i).NtchLn),1),Rotv(RotvRam{i}(Res(i).NtchLn),2),Rotv(RotvRam{i}(Res(i).NtchLn),3),'g','LineWidth',2);
        plot3(Rotv(RotvRam{i}(Res(i).CondLn),1),Rotv(RotvRam{i}(Res(i).CondLn),2),Rotv(RotvRam{i}(Res(i).CondLn),3),'c','LineWidth',2);
    end
    RamCent = mean(Rotv(RotvRam{i},:));
    quiver3(RamCent(1),RamCent(2),RamCent(3),Res(i).RamSymVec(1),Res(i).RamSymVec(2),Res(i).RamSymVec(3),10,'LineWidth',2,'Color','k');
end