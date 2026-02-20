% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [CondLn, RefPtSurf, CondSupVec, CondSymVec, Rmat] = CondLnPrjC3(Rotv, RotvRam, RotfRam, RdgPts, RvN, RCrvns, CrndInd, ChnPt, NchMinInd)
%CondLnPrj returns the line around the upper condylar surface using a method
%of projection of the ridge points associated with the condyle onto various
%planes
%
% Inputs:
%   Rotv
%   RotvRam
%   RotfRam
%   RdgPts
%   RvN
%   RCrvns
%   CrndInd
%   ChnPt
%   NchMinInd
% Outputs:
%   CondLn
%   RefPtSurf
%   CondSupVec
%   CondSymVec
%   Rmat
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%Find condyle ridge seed point - the most distant ridge point from chin
%basal point
RdgPtsI = RdgPts;
%RdgPtsI = cell2mat(RdgPts(:,1));
if isempty(ChnPt)
    YQnt = quantile(Rotv(:,2),3);
    LowThrd = find(Rotv(:,2)<=YQnt(1));
    ChnPt = LowThrd(find(Rotv(LowThrd,3)==max(Rotv(LowThrd,3)),1));
end
RdgChnDst = pdist2(Rotv(RotvRam(RdgPtsI),:),Rotv(ChnPt,:));
RdgSeedPt = RdgPtsI(RdgChnDst==max(RdgChnDst));
RdgSeedPtI = find(RdgChnDst==max(RdgChnDst));

%Find 3D notch minimum and all vertices of the condyle (as defined by notch
%min)
PrjRotv = Rotv(RotvRam,:);
CrndCrd = Rotv(RotvRam(CrndInd),:);
NchMin = Rotv(RotvRam(NchMinInd),:);
%Bmax = max(PrjRotv(PrjRotv(:,3)<NchMin(3),2));
SeedPtDist = pdist2(PrjRotv,PrjRotv(RdgSeedPt,:));
DistThrsh = vecnorm(NchMin-PrjRotv(RdgSeedPt,:),2,2);
InV = find(SeedPtDist<=1.1*DistThrsh);
% DistThrsh = 1.25*DistThrsh - pdist2(PrjRotv(RdgSeedPt,:),mean(PrjRotv(InV,:)));
SeedPtDist = pdist2(PrjRotv,mean(PrjRotv(InV,:)));
InV = find(SeedPtDist<=0.65*DistThrsh); %Parameter change
CondCrds = RvN(InV,:);
% RadVecs = PrjRotv(InV,:)-mean(PrjRotv(InV,:));
% RadVecs = RadVecs./vecnorm(RadVecs,2,2);

%Construct tensor based on unit normal vectors
MatNV = [sum(CondCrds(:,1).^2),sum(CondCrds(:,1).*CondCrds(:,2)),sum(CondCrds(:,1).*CondCrds(:,3));...
    sum(CondCrds(:,2).*CondCrds(:,1)),sum(CondCrds(:,2).^2),sum(CondCrds(:,2).*CondCrds(:,3));...
    sum(CondCrds(:,3).*CondCrds(:,1)),sum(CondCrds(:,3).*CondCrds(:,2)),sum(CondCrds(:,3).^2)];
% MatNR = [sum(RadVecs(:,1).^2),sum(RadVecs(:,1).*RadVecs(:,2)),sum(RadVecs(:,1).*RadVecs(:,3));...
%     sum(RadVecs(:,2).*RadVecs(:,1)),sum(RadVecs(:,2).^2),sum(RadVecs(:,2).*RadVecs(:,3));...
%     sum(RadVecs(:,3).*RadVecs(:,1)),sum(RadVecs(:,3).*RadVecs(:,2)),sum(RadVecs(:,3).^2)];
[eVecs,~] = eig((MatNV));
VecSel = find(abs(eVecs(1,:))==max(abs(eVecs(1,:))));
Rem = setdiff([1,2,3],VecSel);
Rmat = [eVecs(:,VecSel),eVecs(:,Rem(1)),eVecs(:,Rem(2))];
AxDir = [sign(Rmat(1,1)),sign(Rmat(2,2)),sign(Rmat(3,3))];
Rmat = Rmat.*AxDir;



%Use Huersitics to find correct order of eigenvectors
% Ord = perms([1,2,3]);
% for i =1:size(Ord,1)
%     Rmat = [eVecs(:,Ord(i,1)),eVecs(:,Ord(i,2)),eVecs(:,Ord(i,3))];
%     AxDir = [sign(Rmat(1,1)),sign(Rmat(2,2)),sign(Rmat(3,3))];
%     Rmat = Rmat.*AxDir;
%     PrjRotv = (Rmat'* Rotv')';
%     RdgSeedPtCrd = PrjRotv(RotvRam(RdgSeedPt),:);
%     CrndCrd = PrjRotv(RotvRam(CrndInd),:);
%     MaxDims = max(PrjRotv(RotvRam,:))-min(PrjRotv(RotvRam,:));
%     %Test of condyle is "behind" the coranoid, that both are on the "upper" part
%     %and that the the smallest perpendicular dimension is parallel to the x
%     %axis
%     if RdgSeedPtCrd(3) < CrndCrd(3) && RdgSeedPtCrd(2) > 0 && ...
%             CrndCrd(2) > 0 && MaxDims(1) < MaxDims(2) &&...
%             MaxDims(1) < MaxDims(3) && MaxDims(1)/MaxDims(3) < 0.75
%         break
%     end
% end

%Rotate notch minimum and find condylar ridge points, vertices and normals
PrjRotv = (Rmat'* Rotv')';
RdgSeedPtCrd = PrjRotv(RotvRam(RdgSeedPt),:);
CrndCrd = PrjRotv(RotvRam(CrndInd),:);
NchMin = (Rmat'* NchMin')';
YQnt = quantile(PrjRotv(RotvRam,2),100);
if RdgSeedPtCrd(2) < YQnt(65) || CrndCrd(2) < YQnt(65) || NchMin(2) < YQnt(65)
    Rmat = [Rmat(:,1),Rmat(:,3),Rmat(:,2)];
    AxDir = [sign(Rmat(1,1)),sign(Rmat(2,2)),sign(Rmat(3,3))];
    Rmat = Rmat.*AxDir;
    PrjRotv = (Rmat'* Rotv')';
    NchMin = (Rmat'* NchMin')';
end

RCang = RamCondAng(RvN,Rmat);
PrjRvN = (Rmat'* RvN')';
CondRdgPts = RdgPtsI(PrjRotv(RotvRam(RdgPtsI),3)<NchMin(3) &...
    PrjRotv(RotvRam(RdgPtsI),2)>=NchMin(2));
% InV = find(PrjRotv(RotvRam,3) < NchMin(3) & PrjRotv(RotvRam,2) > NchMin(2));

%reindex subset of triangles
RotvCndInd = InV;
RotfCndInd = find(any(ismember(RotfRam,RotvCndInd),2));
RotfCnd = RotfRam(RotfCndInd,:);
[InV,~,IndUptN] = unique(RotfCnd);
RotfCnd = reshape(IndUptN,size(RotfCnd));


%Initial reference point X&Z average on plane; Y top projection of XZ
RefPt = [mean(PrjRotv(RotvRam(InV),1)),mean(PrjRotv(RotvRam(InV),3))];
RotfCndCrd = zeros(3,3,size(RotfCnd,1));
for i=1:size(RotfCnd,1)
    RotfCndCrd(:,:,i) = PrjRotv(RotvRam(InV(RotfCnd(i,:))),:);
end
YTris = find(squeeze(~all(RotfCndCrd(:,1,:)<=RefPt(1),1) &...
    ~all(RotfCndCrd(:,1,:)>=RefPt(1),1) &...
    ~all(RotfCndCrd(:,3,:)<=RefPt(2),1) &...
    ~all(RotfCndCrd(:,3,:)>=RefPt(2),1)));
YVals = squeeze(RotfCndCrd(:,2,YTris));
% BHght = mean(YVals(YVals>mean(PrjRotv(RotvRam(InV),2))),'all');
BHght = max(mean(YVals));
RefPt = [RefPt(1),BHght,RefPt(2)];
RefPtSurf = RefPt;

% theta = atan2(PrjRotv(RotvRam(CondRdgPts),3)-RefPt(3),PrjRotv(RotvRam(CondRdgPts),1)-RefPt(1));
% PointsXYZ = [theta,PrjRotv(RotvRam(CondRdgPts),:)];

%calculates a single "shadow point" along a closed circle (projected on the plane) 
%Finds the shadow point in 4 degrees intevals (slices) using 5 different
%height for the reference point. For more details on finding the point see
%FindCondEdgPts (finding the tangent) and EdgPtCent (centering to a single
%point)
CentPts = zeros(90,4,15);
i = 1;
HFac = 1;
EdgePts = [];
while HFac <= 3*(BHght-mean(PrjRotv(RotvRam(InV),2)))
    HFac = 1+ i*0.25;  
    RefPt(2) = BHght+HFac;
    EdgePts = FindCondEdgPts(RefPt,PrjRotv,RotvRam,InV,PrjRvN,RotfCnd);
    % EdgePts = [EdgePts;EdgePtsi];
    CentPts(:,:,i) = EdgPtCent(RefPt,PrjRotv,RotvRam,EdgePts);
    i=i+1;
end


%Remove noisy results heights. If reference point is too low, given that there is no
%procedure for removing points blocked by the surface, results can be
%extremely messy.
CentDists = squeeze(vecnorm(CentPts(:,2:4,:) - RefPt,2,2));
RefHFltr = ~isoutlier(std(CentDists),'quartiles');
CentPts = CentPts(:,:,RefHFltr);
CentDists = squeeze(vecnorm(CentPts(:,2:4,:) - RefPt,2,2));

%Another round of cleaning
CentHVar = std(CentDists,[],2);
Otlrs = find(isoutlier(CentHVar));
Inlrs = setdiff((1:90)',Otlrs);
CentPtsM = zeros(90,3);
CentPtsM(Inlrs,:) = mean(CentPts(Inlrs,2:4,:),3);
DistM = vecnorm(CentPtsM-RefPt,2,2);
DistM = [DistM;DistM(1)];
DistM = filloutliers(DistM,"pchip");
for i=1:size(Otlrs,1)
    CentPtsM(Otlrs(i),:) = CentPts(Otlrs(i),2:4,...
        find(abs(CentDists(Otlrs(i),:)-DistM(Otlrs(i)))==min(abs(CentDists(Otlrs(i),:)-DistM(Otlrs(i)))),1));
end
CondLn = [CentPtsM;CentPtsM(1,:)];


CondLn = (Rmat * CondLn')';
RefPtSurf = (Rmat * RefPtSurf')';
CondSupVec = (Rmat * [0,1,0]')';
CondSymVec = (Rmat * [1,0,0]')';
%CondAng = acos(dot(CondVec,[0,1,0],2));


% Cxyz = fmincon(@(cent) RefPtOptFnc(cent,PrjRotv,RotvRam,InV,PrjRvN,RotfCnd,CondRdgPts),...
% RefPt,[],[],[],[],[min(PointsXYZ(:,2)),BHght+0.1*(BHght-min(PrjRotv(RotvRam(InV),2))),min(PointsXYZ(:,4))],...
% [max(PointsXYZ(:,2)),BHght+2*(BHght-min(PrjRotv(RotvRam(InV),2))),max(PointsXYZ(:,4))],...
% @(x)ConsT(x,RefPt,PointsXYZ),optimoptions("fmincon","Algorithm","interior-point",...
%     "SubproblemAlgorithm","cg"));
% 





% %Get Curvness values for condylar ridge points, calculate the centroid and
% %find theta
% RCrvTh = quantile(RCrvns(RdgPtsI(CondRdgPts2)),100);
% CondRdgPtsCrv = CondRdgPts2(RCrvns(RdgPtsI(CondRdgPts2))>RCrvTh(1));
% RdgPtsCrv = RCrvns(RdgPtsI(CondRdgPtsCrv));
% centroid = mean(PrjRotv(RotvRam(RdgPtsI(CondRdgPtsCrv)),:));
% %centroid = RdgSeedPt;
% 
% theta = atan2(PrjRotv(RotvRam(RdgPtsI(CondRdgPtsCrv)),3)-centroid(3),PrjRotv(RotvRam(RdgPtsI(CondRdgPtsCrv)),1)-centroid(1));
% 
% %move to origin and convert to spherical coordinates. PointsTP contains the
% %independent variable angles theta (0 -> pi) and phi (-pi -> pi) while
% %PointsR contains the dependent radius
% PointsXYZ = [theta,PrjRotv(RotvRam(RdgPtsI(CondRdgPtsCrv)),:)];
% PointsO = PointsXYZ(:,2:4) - centroid;
% PointsR = vecnorm(PointsO,2,2);
% PointsTP = acos(PointsO(:,2)./PointsR);
% PointsTP = [PointsTP, atan2(PointsO(:,3),PointsO(:,1))];
% 
% %Assigns weights to the condyle ridge points based on their postion - with
% %preference towards posterior and superior points
% PosFac = abs(PointsXYZ(:,3).*PointsXYZ(:,4));
% Weight = PosFac./sum(PosFac);
% 
% %Optimizes the position of the centroid by maximizing the wieghted standard
% %deviation of phi (with constarins)
% 
% Cxyz = fmincon(@(cent) StdPhiBins(cent,PointsXYZ,Weight),...
% centroid,[],[],[],[],[min(PointsXYZ(:,2)),0,min(PointsXYZ(:,4))],[max(PointsXYZ(:,2)),100,max(PointsXYZ(:,4))],...
% @(x)ConsT(x,centroid,PointsXYZ),optimoptions("fmincon",'Display','none'));
% 
% % Cy = fmincon(@(cent) 1./std(acos((PointsXYZ(:,3)-cent(2))./vecnorm(PointsXYZ(:,2:4)-cent,2,2))),...
% % centroid,[],[],[],[],[min(PointsXYZ(:,2)),0,min(PointsXYZ(:,4))],[max(PointsXYZ(:,2)),100,max(PointsXYZ(:,4))],...
% % @(x)ConsT(x,centroid,PointsXYZ),optimoptions("fmincon",'Display','none'));
% 
% PointsO = PointsXYZ(:,2:4) - Cxyz;
% PointsR = vecnorm(PointsO,2,2);
% PointsTP = acos(PointsO(:,2)./PointsR);
% PointsTP = [PointsTP, atan2(PointsO(:,3),PointsO(:,1))];
% 
% %fit a harmonic function to find a relation between phi and theta. Uses
% %ransac to ignore noise.
% 
% Weight = 0.2.*RdgPtsCrv./sum(RdgPtsCrv)+0.8*Weight;
% [TPModelFun, ModelCons, InlrsCons] = RanSac_Harmonic7(500, 20, 0.05, PointsTP, Weight);
% 
% %fits a spherical harmonics of l=3 (20 coefficients including beta m=0)
% %than runs ransac on that fit onto the selected points
% Weight = 0.2.*RdgPtsCrv(InlrsCons)./sum(RdgPtsCrv(InlrsCons))+0.8*abs(PointsXYZ(InlrsCons,3).*PointsXYZ(InlrsCons,4))./sum(abs(PointsXYZ(InlrsCons,3).*PointsXYZ(InlrsCons,4)));
% [ModelFix, ShprInlrsCons] = RanSac_SphrHar(125, 30, 0.4, [PointsTP(InlrsCons,:),PointsR(InlrsCons)], Weight);
% 
% % SphDistFun = @(ModelShr) abs(SphrHar(ModelShr,PointsTP(:,1:2))-PointsR).*Weight;
% % ModelFix = lsqnonlin(SphDistFun,ones(10,2),repmat(-10,10,2),repmat(10,10,2),optimoptions("lsqnonlin",'Display','none',"UseParallel",true));
% % 
% %  'fitting Spherical Harmonics'
% % toc
% 
% 
% 
% % PointsTPSel = PointsTP(InlrsCons,:);
% % PointsRSel = PointsR(InlrsCons,:);
% % WeightSel = Weight(InlrsCons);
% % 
% 
% % the line in specrical coordinates [theta as predicted by the TPfun, phi
% % -pi -> pi, and R as prediction of the spherical harmonics for the given
% % theta and phi]
% TPR = [TPModelFun(ModelCons,(-pi:2*pi/99:pi)'),(-pi:2*pi/99:pi)',...
%     SphrHar(ModelFix,[TPModelFun(ModelCons,(-pi:2*pi/99:pi)'),(-pi:2*pi/99:pi)'])];
% 
% %convert back to Eulidean coordinates and translate back to place
% X = TPR(:,3).*sin(TPR(:,1)).*cos(TPR(:,2));
% Y = TPR(:,3).*cos(TPR(:,1));
% Z = TPR(:,3).*sin(TPR(:,1)).*sin(TPR(:,2));
% XYZ = [X,Y,Z];
% CondLn = XYZ+Cxyz;
% 
% 
% %Rotate back to mandible's table plane
% % Rmat = RMatAxAng(Rvec,Rang);
% CondLn = (Rmat* (CondLn-centroid)')'+centroid;
end