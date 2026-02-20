% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [NchLn] = NotchLine2(Rotv, RotvRam, CrndInd, NchMinInd, RCrvPts, RCrvns, CondLn, RKDmin, RvN)
%
% Inputs:
%   Rotv
%   RotvRam
%   CrndInd
%   NchMinInd
%   RCrvPts
%   RCrvns
%   CondLn
%   RKDmin
%   RvN
% Outputs:
%   NchLn
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

    %Leave only ridge points that are in the general (+-1/4*pi) direction of 
    %the notch minimum
    %First for the porjection onto the coronoid plane
    CrndvN = RvN(CrndInd,:);
    RdgPtsInd = RCrvPts;
    RdgPtsMov = {RdgPtsInd};
    NchMinCrd = Rotv(RotvRam(NchMinInd),:);
    NchMinVec = (NchMinCrd - Rotv(RotvRam(CrndInd),:));
    NchMinDst = dot(NchMinVec,CrndvN,2);
    RdgPtsVec = (Rotv(RotvRam(RdgPtsInd),:) - Rotv(RotvRam(CrndInd),:));
    RdgPtsDst = dot(RdgPtsVec,repmat(CrndvN,size(RdgPtsVec,1),1),2);
    NchMinPrj = NchMinCrd - NchMinDst*CrndvN;
    RdgPtsPrj = Rotv(RotvRam(RdgPtsInd),:) - RdgPtsDst*CrndvN;
    NchMinPNVec = (NchMinPrj - Rotv(RotvRam(CrndInd),:))./vecnorm((NchMinPrj - Rotv(RotvRam(CrndInd),:)),2,2);
    RdgPtsPNVec = (RdgPtsPrj - Rotv(RotvRam(CrndInd),:))./vecnorm((RdgPtsPrj - Rotv(RotvRam(CrndInd),:)),2,2);
    RdgPtsPrjDir = acos(dot(repmat(NchMinPNVec,size(RdgPtsPNVec,1),1),RdgPtsPNVec,2));

    %Next in 3D general direction
    NchMinNVec = NchMinVec./vecnorm(NchMinVec,2,2);
    RdgPtsNVec = RdgPtsVec./vecnorm(RdgPtsVec,2,2);
    RdgPtsDir = acos(dot(repmat(NchMinNVec,size(RdgPtsNVec,1),1),RdgPtsNVec,2));
   
    RdgPtsInd = RdgPtsInd(any([RdgPtsPrjDir<=0.25*pi,RdgPtsDir<=0.25*pi],2));    %Factor! Ridge point direction filter. default = 0.25*pi
    
    %Fit a line direction to the ridge points closest to the start point. Start
    %with 0.15 distance and increase until there are at least 15 points to fit
    CorCondDist = pdist2(Rotv(RotvRam(CrndInd),:),Rotv(RotvRam(CondLn),:));
    TrgtPt = CondLn(find(CorCondDist == min(CorCondDist),1));
    TrgtDir = Rotv(RotvRam(TrgtPt),:) - Rotv(RotvRam(CrndInd),:);
    TrgtDir = TrgtDir./ vecnorm(TrgtDir,2,2);
    OTrgtDir = TrgtDir;
    RdgPtsDist = pdist2(Rotv(RotvRam(RdgPtsInd),:),Rotv(RotvRam(CrndInd),:));
    % d = 0.15;
    % FitSmp = RdgPtsInd(RdgPtsDist<=d);
    % while numel(FitSmp)<10
    %     d=d+0.05;
    %     FitSmp = RdgPtsInd(RdgPtsDist<=d);
    %     FitSmpDir = Rotv(RotvRam(FitSmp),:) - Rotv(RotvRam(CrndInd),:);
    %     FitSmpDir = FitSmpDir./vecnorm(FitSmpDir,2,2);
    %     FitSmpAng = acos(dot(repmat(TrgtDir,size(FitSmpDir,1),1),FitSmpDir,2));
    %     FitSmp = FitSmp(FitSmpAng<0.3);
    % end
    % RdgPtsVec = (Rotv(RotvRam(FitSmp),:) - Rotv(RotvRam(CrndInd),:))...
    %     ./vecnorm(Rotv(RotvRam(FitSmp),:) - Rotv(RotvRam(CrndInd),:),2,2);
    % %Find leat square vector fit (1st principal component), rotate if necessary
    % PcDirVec = pca(Rotv(RotvRam(FitSmp),:),'Weights',RCrvns(FitSmp)./sum(RCrvns(FitSmp)));
    % PcDirVec = PcDirVec(:,1)';
    % FitDirVec = mean(RdgPtsVec)./vecnorm(mean(RdgPtsVec),2,2);
    % if acos(dot(FitDirVec,PcDirVec,2)) >0.5*pi
    %     PcDirVec = -PcDirVec;
    % end
    
    %Second point in line is that direction distance threshold away
    DirPt = Rotv(RotvRam(CrndInd),:) + 0.25*TrgtDir;
    PcDirVecPrv = TrgtDir;
    RdgPtsInd = setdiff(RdgPtsInd,RdgPtsInd(RdgPtsDist<=0.25));
    NchLn = [CrndInd,RdgPtsInd(knnsearch(Rotv(RotvRam(RdgPtsInd),:),DirPt,'K',1))];
    RdgPtsInd = setdiff(RdgPtsInd,NchLn);
    RdgPtsMov = [RdgPtsMov;{RdgPtsInd}];
    TDist = pdist2(Rotv(RotvRam(CrndInd),:),Rotv(RotvRam(CondLn),:));
    TDist = min(TDist);
    CorCondDist = pdist2(Rotv(RotvRam(NchLn(end)),:),Rotv(RotvRam(CondLn),:));
    TrgtPt = CondLn(find(CorCondDist == min(CorCondDist),1));
    TrgtDir = Rotv(RotvRam(TrgtPt),:) - Rotv(RotvRam(NchLn(end)),:);
    TrgtDir = TrgtDir./ vecnorm(TrgtDir,2,2);
    RelProg = min(CorCondDist)./TDist;
    CorCondDistOld = min(CorCondDist);
   
    %start loop
    LoopFlg = true;
    
    
    while LoopFlg
        RdgPtsDist = pdist2(Rotv(RotvRam(RdgPtsInd),:),Rotv(RotvRam(NchLn(end)),:));
        d = 0.15;
        FitSmp = RdgPtsInd(RdgPtsDist<=d);
        %Gradually increase distance if not enough points for fitting but only
        %up to a threshold, otherwise end loop. The threshold has a factor
        %derived from the progress along the line, as it approaches the end,
        %less distance is allowed
        while numel(FitSmp)<15 && d <=0.50 + 4.5*RelProg
            d=d+0.05;
            FitSmp = RdgPtsInd(RdgPtsDist<=d);
            %Clean FitSmp points "behind" the line end-point
            % RdgPtsVec = (Rotv(RotvRam(FitSmp),:) - Rotv(RotvRam(NchLn(end)),:))...
            %     ./vecnorm(Rotv(RotvRam(FitSmp),:) - Rotv(RotvRam(NchLn(end)),:),2,2);
            % RdgPtsPrjDir = dot(repmat(TrgtDir,size(RdgPtsVec,1),1),RdgPtsVec,2);
            % RdgPtsInd = setdiff(RdgPtsInd,FitSmp(RdgPtsPrjDir<0));
            % FitSmp = setdiff(FitSmp,FitSmp(RdgPtsPrjDir<0));
            % RdgPtsDist = pdist2(Rotv(RotvRam(RdgPtsInd),:),Rotv(RotvRam(NchLn(end)),:));
            
        end 
        if numel(FitSmp)<10
            break
        end
    %     %Calculate fit points directions and the line (PC) direction, check
    %     %which is best in line with previous direction, rotate if needed

        %Curvedness used as weights
        PcDirVec = pca(Rotv(RotvRam(FitSmp),:),'Weights',RCrvns(FitSmp)./sum(RCrvns(FitSmp)));
%         RdgDirVec = RKDmin(FitSmp,:);
%         Dist = dot(repmat(PcDirVecPrv,size(RdgDirVec,1),1),RdgDirVec,2);
%         RdgDirVec = mean(Dist.*RdgDirVec)./vecnorm(mean(Dist.*RdgDirVec),2,2);
        PcDirVecAngs = acos(dot(PcDirVec',repmat(TrgtDir,3,1),2));
        CorVec = find(abs(PcDirVecAngs-0.5*pi)==max(abs(PcDirVecAngs-0.5*pi)));
        PcDirVec = PcDirVec(:,CorVec)';
    %     FitDirVec = mean(RdgPtsVec)./vecnorm(mean(RdgPtsVec),2,2);
        %reverse if opposite to the direction of previous section
        if acos(dot(TrgtDir,PcDirVec,2)) >0.5*pi
            PcDirVec = -PcDirVec;
        end
        % include ridge direction (reflected by direction of min curvature)
        RdgDirVec = RKDmin(FitSmp,:);
        Flp = acos(dot(repmat(TrgtDir,size(RdgDirVec,1),1),RdgDirVec,2)) >0.5*pi;
        RdgDirVec(Flp,:) = -RdgDirVec(Flp,:);
        Dist = dot(repmat(PcDirVec,size(RdgDirVec,1),1),RdgDirVec,2);
        RdgDirVec = mean(Dist.*RdgDirVec)./vecnorm(mean(Dist.*RdgDirVec),2,2);
        PcDirVec = sum([0.2.*PcDirVec;0.15.*RdgDirVec;0.3.*PcDirVecPrv;0.25.*TrgtDir;0.1.*OTrgtDir]); %[0.2,0.25,0.2,0.15,0.2],
        PcDirVec = PcDirVec./vecnorm(PcDirVec,2,2);
    
        %Find the next point using the vector, snap to nearest ridge point if
        %under a distance threshold, otherwise break
        DirPt = Rotv(RotvRam(NchLn(end)),:)+0.25*PcDirVec;
    
        RdgPtDist = pdist2(Rotv(RotvRam(FitSmp),:),DirPt);
%         RdgPtQal = RdgPtDist <=0.5;
%         if any(RdgPtQal)
        NchLn = [NchLn,FitSmp(RdgPtDist==min(RdgPtDist))];
        if numel(NchLn)<10
            NchLnEndSeg = NchLn(2:end);
        else
            NchLnEndSeg = NchLn((end-8):(end));
        end
        PcDirVecPrv = pca(Rotv(RotvRam(NchLnEndSeg),:),'Weights',RCrvns(NchLnEndSeg)./sum(RCrvns(NchLnEndSeg)));
        PcDirVecPrv = PcDirVecPrv(:,1)';
        if dot(PcDirVecPrv,Rotv(RotvRam(NchLn(end)),:)-Rotv(RotvRam(NchLn(end-1)),:),2) <0
            PcDirVecPrv = -PcDirVecPrv;
        end
        
        RdgPtsInd = setdiff(RdgPtsInd,RdgPtsDist<=0.25);
        RdgPtsInd = setdiff(RdgPtsInd,NchLn);
        RdgPtsDist = pdist2(Rotv(RotvRam(RdgPtsInd),:),Rotv(RotvRam(NchLn(end)),:));

        % RdgPtsVec = (Rotv(RotvRam(RdgPtsInd),:) - Rotv(RotvRam(NchLn(end)),:))...
        %     ./vecnorm(Rotv(RotvRam(RdgPtsInd),:) - Rotv(RotvRam(NchLn(end)),:),2,2);
        RdgPtsVec = (Rotv(RotvRam(RdgPtsInd),:) - mean(Rotv(RotvRam(NchLnEndSeg),:)));
        RdgPtsPrjDir = dot(repmat(PcDirVecPrv,size(RdgPtsVec,1),1),RdgPtsVec,2);
        RdgPtsInd = setdiff(RdgPtsInd,RdgPtsInd(RdgPtsPrjDir>=-0.5 & RdgPtsPrjDir<=0.5 & RdgPtsDist <= d));
        RdgPtsMov = [RdgPtsMov;{RdgPtsInd}];
%         else
%             LoopFlg = false;
%         end
        
        
        %If current end point is less than threshold distance to the condylar
        %line, stop adding and snap to it
        CorCondDist = pdist2(Rotv(RotvRam(NchLn(end)),:),Rotv(RotvRam(CondLn),:));
        RelProg = min(CorCondDist)./TDist;
        TrgtPt = CondLn(find(CorCondDist == min(CorCondDist),1));
        TrgtDir = Rotv(RotvRam(TrgtPt),:) - Rotv(RotvRam(NchLn(end)),:);
        TrgtDir = TrgtDir./ vecnorm(TrgtDir,2,2);
        if min(CorCondDist)<=0.5
            NchLn = [NchLn,TrgtPt];
            LoopFlg = false;
        elseif RelProg <= 0.25 && min(CorCondDist) > CorCondDistOld 
            LoopFlg = false;
        else
            CorCondDistOld = min(CorCondDist);
        end
    end

end
%%
% %Find ridge points in cylinder between the two points. calculate distance
%     %from the main line
%     InCylPts = PtsInCyl(Rotv(RotvRam(NchLn(end-1)),:),Rotv(RotvRam(NchLn(end)),:),Rotv(RotvRam(RdgPtsInd),:),1.5);
%     if numel(InCylPts)>1
%         InCylDst = vecnorm(cross(repmat(Rotv(RotvRam(NchLn(end)),:)-Rotv(RotvRam(NchLn(end-1)),:),[size(InCylPts,1),1]),(Rotv(RotvRam(NchLn(end-1)),:)-Rotv(RotvRam(RdgPtsInd(InCylPts)),:))),2,2)...
%             /vecnorm(Rotv(RotvRam(NchLn(end)),:)-Rotv(RotvRam(NchLn(end-1)),:),2,2);
%         RdgPtsInd = setdiff(RdgPtsInd,InCylPts);
%     
%         %remove cylinder end point (also included in list, distance = 0), get these
%         %points ridge directions, reverse if neeeded, caluclate weighted sum
%     %     InCylPts = InCylPts(InCylDst~=0);
%     %     InCylDst = InCylDst(InCylDst~=0);
%         InCylRdgDir = KDmin(InCylPts,:);
%         Rev = acos(dot(repmat(Rotv(RotvRam(NchLn(end)),:)-Rotv(RotvRam(NchLn(end-1)),:),[size(InCylPts,1),1]),InCylRdgDir,2))>0.5*pi;
%         InCylRdgDir(Rev,:) = -InCylRdgDir(Rev,:);
%         InCylWght = abs(1.5-InCylDst)./sum(abs(1.5-InCylDst));
%         CylDirVec = sum(InCylWght.*InCylRdgDir)./vecnorm(sum(InCylWght.*InCylRdgDir),2,2);
%     elseif numel(InCylPts) == 1
%         CylDir = (Rotv(RotvRam(NchLn(end)),:)-Rotv(RotvRam(NchLn(end-1)),:));
%         CylDir = CylDir./vecnorm(CylDir,2,2);
%         CylDirVec = mean([KDmin(InCylPts,:);CylDir])./vecnorm(mean([KDmin(InCylPts,:);CylDir]),2,2);
%     else
%         CylDir = (Rotv(RotvRam(NchLn(end)),:)-Rotv(RotvRam(NchLn(end-1)),:));
%         CylDirVec = CylDir./vecnorm(CylDir,2,2);
%     end
%     
%     %average with end point ridge direction, reverse if needed
%     if acos(dot(KDmin(NchLn(end),:),CylDirVec,2)) > 0.5*pi
%         LnDirVec = mean([CylDirVec;-KDmin(NchLn(end),:)])...
%             ./vecnorm(mean([CylDirVec;-KDmin(NchLn(end),:)]),2,2);
%     else
%         LnDirVec = mean([CylDirVec;KDmin(NchLn(end),:)])...
%             ./vecnorm(mean([CylDirVec;KDmin(NchLn(end),:)]),2,2);
%     end
%     DirPt = Rotv(RotvRam(NchLn(end)),:)+0.5*LnDirVec;
%     
%     %neighboring triangles to current line end amd their normals
%     PotTriPrj = vertexAttachments(TrR,NchLn(end));
%     PotTriPrj = unique(neighbors(TrR,unique(neighbors(TrR,cell2mat(PotTriPrj)'))));
%     PotTriNV = faceNormal(TrR,PotTriPrj);
%     PotTriCC = circumcenter(TrR,PotTriPrj);
%     PrjDirVecs = DirPt-PotTriCC;
%     PrjDists = dot(PrjDirVecs,PotTriNV,2);
%     DirPtPrj = PotTriCC - PrjDists.*PotTriNV;
% 
%     RdgPtDist = pdist2(Rotv(RotvRam(RdgPtsInd),:),mean(DirPtPrj));
%     RdgPtQal = RdgPtDist <=0.5;
%     if any(RdgPtQal)
%         NchLn = [NchLn,RdgPtsInd(RdgPtDist==min(RdgPtDist))];
%         RdgPtsInd = setdiff(RdgPtsInd,NchLn);
%     else
%         LoopFlg = false;
%     end
%     CorCondDist = pdist2(Rotv(RotvRam(NchLn(end)),:),Rotv(RotvRam(CondLn),:));
%     
%     if min(CorCondDist)<=0.5
%         TrgtPt = CondLn(find(CorCondDist == min(CorCondDist),1));
%         NchLn = [NchLn,CondLn(TrgtPt)];
%         LoopFlg = false;
%     end

% %Start line in direction of ridge from the coronoid. Check if not reverse
% %to direction to notch minimum, and if so reverse. Advance 0.5mm and snap
% %to nearest ridge point
% 
% if acos(dot(NchMinVec,RKDmin(CrndInd,:),2)) > 0.5*pi
%     LnDirVec = -RKDmin(CrndInd,:);
% else
%     LnDirVec = RKDmin(CrndInd,:);
% end
% DirPt = Rotv(RotvRam(CrndInd),:)+0.5*LnDirVec;
% NchLn = [CrndInd,RdgPtsInd(knnsearch(Rotv(RotvRam(RdgPtsInd),:),DirPt,'K',1))];
% RdgPtsInd = setdiff(RdgPtsInd,NchLn);
% LoopFlg = true;


% %Find nearest neigboring curvature points culter whose planar projection
% %angle (on Coronoid plane) is closest to the mean between the ridge
% %direction and the vector pointing towards the nothc minimum point
% 
% CrndRdgNN = knnsearch(Rotv(RotvRam(RdgPtsInd),:),Rotv(RotvRam(CrndInd),:),'K',200);
% CrndRdgNN = RdgPtsInd(CrndRdgNN);
% NNv =  Rotv(RotvRam(CrndRdgNN),:) - Rotv(RotvRam(CrndInd),:); 
% NNdists = dot(NNv,repmat(CrndvN,size(NNv,1),1),2);
% PPrj = Rotv(RotvRam(CrndRdgNN),:) - NNdists*CrndvN;
% theta = atan2(PPrj(:,3)-Rotv(RotvRam(CrndInd),3),PPrj(:,1)-Rotv(RotvRam(CrndInd),1));
% Lnks = linkage(theta,"ward");
% DistDiffs = (Lnks(2:end,3)-Lnks(1:end-1,3))./Lnks(1:end-1,3);
% CutDist = Lnks(find(DistDiffs(end-3)==max(DistDiffs(end-3)))+(size(DistDiffs,1)-2)+1,3);
% CutDist = CutDist - CutDist.*0.05;
% Clsts = cluster(Lnks,"criterion","distance","cutoff",CutDist);
% 
% %direction vectors of clusters
% ClstVecs = zeros(numel(unique(Clsts)),3);
% for j=1:numel(unique(Clsts))
%     ClstVecs(j,:) = mean(PPrj(Clsts==j,:)-Rotv(RotvRam(CrndInd),:))./...
%         vecnorm(mean(PPrj(Clsts==j,:)-Rotv(RotvRam(CrndInd),:)),2,2); 
% end
