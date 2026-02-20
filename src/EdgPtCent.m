% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [CentPt] = EdgPtCent(RefPt,PrjRotv,RotvRam,EdgePts)
    %Divide into bins
%
% Inputs:
%   RefPt
%   PrjRotv
%   RotvRam
%   EdgePts
% Outputs:
%   CentPt
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.
    PhiInt = linspace(-pi,pi,91);
    PointsO = PrjRotv(RotvRam(EdgePts),:) - RefPt;
    EdgeThet = atan2(PointsO(:,3),PointsO(:,1));
    ThetBin = discretize(EdgeThet,PhiInt);
    
    
    %Calculate colsest
    CentPt = zeros(90,4);
    for i=1:90
        InSlcInd = find((ThetBin==i));
        InSlcDists = vecnorm(PointsO(ThetBin==i,:),2,2);
        Clst = InSlcInd(InSlcDists<=min(InSlcDists)+1);
        CentPt(i,:) = [mean([PhiInt(i);PhiInt(i+1)]),mean(PrjRotv(RotvRam(EdgePts(Clst)),:),1)];
    end
    
    %Clean results
    % for i=1:5
    %     %Local
    %     PtDist = vecnorm(CentPt(:,2:4)-RefPt,2,2);
    %     OLrs = find(isoutlier(PtDist,'movmedian',5));
    %     OLVecs = (CentPt(OLrs,2:4)-RefPt)./vecnorm(CentPt(OLrs,2:4)-RefPt,2,2);
    %     PtDist = filloutliers(PtDist,"pchip","movmedian",5);
    %     CentPt(OLrs,2:4) = RefPt + OLVecs.*PtDist(OLrs);
    %     %Global
    %     OLrs = find(isoutlier(PtDist));
    %     OLVecs = (CentPt(OLrs,2:4)-RefPt)./vecnorm(CentPt(OLrs,2:4)-RefPt,2,2);
    %     PtDist = filloutliers(PtDist,"pchip");
    %     CentPt(OLrs,2:4) = RefPt + OLVecs.*PtDist(OLrs);
    % end
    PtDist = vecnorm(CentPt(:,2:4)-RefPt,2,2);
    OLrs = find(isoutlier(PtDist,'movmedian',5));
    OLrs = [OLrs;find(isoutlier(PtDist))];
    CentPt(OLrs,2:4) = NaN;
    CentPt = fillmissing(CentPt,"pchip",1);
end