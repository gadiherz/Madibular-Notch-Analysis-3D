% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [NchMinInd,CrndInd] = CrndNchMin(Rotv,RotvRam,RCrvns,RKGauss)
%CrndNchMin Finds the index of the coronoid and notch minimum
%   it uses the iterative peak curvature detection procedure to detect the
%   2D coords (Z and Y) of the points, then looks for the right point on
%   mesh using the Curvdness and Gauss curvature values of nearby points
%
% Inputs:
%   Rotv
%   RotvRam
%   RCrvns
%   RKGauss
% Outputs:
%   NchMinInd
%   CrndInd
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%Find boundary about projection onto the symmetry plane
b4Pos = boundary(Rotv(RotvRam,3),Rotv(RotvRam,2));
b4Pos = [Rotv(RotvRam(b4Pos),3),Rotv(RotvRam(b4Pos),2)];

% Strong smoothing dfor curvature, weak for coordinates
[NewDist,~,~,~,~,K,~,~,~] = FourierCoords(b4Pos,6);
[~,~,~,~,~,~,Xpred,Ypred,~] = FourierCoords(b4Pos,100);

%Find minima points
PeaksNeg = find(all([islocalmin(K),sign(K)<0],2));
NchMinZY = [Xpred(PeaksNeg(Ypred(PeaksNeg)==max(Ypred(PeaksNeg)))),Ypred(PeaksNeg(Ypred(PeaksNeg)==max(Ypred(PeaksNeg))))];


%Find maxima points
PeaksPos = find(all([islocalmax(K),sign(K)>0,Xpred>NchMinZY(1),Ypred>NchMinZY(2)],2));
PeaksPos = PeaksPos(K(PeaksPos)==max(K(PeaksPos)));
PeaksNeg = PeaksNeg(Ypred(PeaksNeg)==max(Ypred(PeaksNeg)));

%Refines the location of points based on lower smoothing
[NewDist,~,~,~,~,K,~,~,~] = FourierCoords(b4Pos,100);

PeaksNeg = find(all([islocalmin(K),sign(K)<0,[zeros(PeaksNeg-20,1);ones(40,1);zeros(numel(NewDist)-(PeaksNeg-20+40),1)] ],2));
PeaksNeg = PeaksNeg(K(PeaksNeg)==min(K(PeaksNeg)));
NchMinZY = [Xpred(PeaksNeg), Ypred(PeaksNeg)];
PeaksPos = find(all([islocalmax(K),sign(K)>0,Xpred>NchMinZY(1),Ypred>NchMinZY(2),[zeros(PeaksPos-20,1);ones(40,1);zeros(numel(NewDist)-(PeaksPos-20+40),1)]],2));
CrndZY = [Xpred(PeaksPos(Ypred(PeaksPos)==max(Ypred(PeaksPos)))),Ypred(PeaksPos(Ypred(PeaksPos)==max(Ypred(PeaksPos))))];

%Snap to 3D with highest curvedness 
PtDists = pdist2(NchMinZY,[Rotv(RotvRam,3),Rotv(RotvRam,2)]);
Snp = find(PtDists == min(PtDists),1);
NchMinZY = [Rotv(RotvRam(Snp),3),Rotv(RotvRam(Snp),2)];
PotNchMin = find(Rotv(RotvRam,2)>=NchMinZY(2)-0.25 & Rotv(RotvRam,2)<=NchMinZY(2)+0.25...
    & Rotv(RotvRam,3)>=NchMinZY(1)-0.25 & Rotv(RotvRam,3)<=NchMinZY(1)+0.25);
NchMinInd = PotNchMin(filloutliers(RCrvns(PotNchMin),'center')==max(filloutliers(RCrvns(PotNchMin),'center')));
if size(NchMinInd,1)>1
    PtDists = pdist2(mean(Rotv(RotvRam(NchMinInd),:)),Rotv(RotvRam,:));
    Snp = find(PtDists == min(PtDists),1);
    NchMinInd = Snp;
end

%Snap to 3D with highest Gauss
PtDists = pdist2(CrndZY,[Rotv(RotvRam,3),Rotv(RotvRam,2)]);
Snp = find(PtDists == min(PtDists),1);
CrndZY = Rotv(RotvRam(Snp),:);
PotCrnd = find(Rotv(RotvRam,2)>=CrndZY(2)-0.5...% & Rotv(RotvRam,2)<=CrndZY(2)+0.4...
    & Rotv(RotvRam,3)>=CrndZY(3)-0.5 & Rotv(RotvRam,3)<=CrndZY(3)+0.5);
SnpDist = pdist2(CrndZY,Rotv(RotvRam(PotCrnd),:));
SnpPtFix = mink(SnpDist,2);
SnpDist(SnpDist==0) = SnpPtFix(2);
SnpDist = 1./SnpDist; %inverse distance from snap point
SnpDist = (SnpDist./max(SnpDist))'; %Normalized to 0-1
SnpDir = (Rotv(RotvRam(PotCrnd),1) - CrndZY(1)).*sign(CrndZY(1)); %Direction of difference, towards the buccal is positive 
SnpDir = SnpDir - min(SnpDir);
SnpDir = SnpDir ./ max(SnpDir); %Normalized to 0-1
SnpGauss = filloutliers(RKGauss(PotCrnd),'center');
SnpGauss = SnpGauss - min(SnpGauss);
SnpGauss = SnpGauss ./ max(SnpGauss); %Normalized to 0-1
CstFun = 0.5.*SnpGauss + 0.3.*SnpDir + 0.2.*SnpDist;

CrndInd = PotCrnd(CstFun==max(CstFun)); 

end