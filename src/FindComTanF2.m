% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Ang, PtInd] = FindComTanF(Peaks,NewDist, Theta, K,Xpred, Ypred)
%FindComTanF Finds the points through which the best common tangent (not
%intersecting the curve) passes for a section. 
%   It looks for points whos tanget is the same line (having the same angle
%   and inersection with Y at X=0).
%
% Inputs:
%   Peaks
%   NewDist
%   Theta
%   K
%   Xpred
%   Ypred
% Outputs:
%   Ang
%   PtInd
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%Recieves the two peaks and finds the section indices by order - important
%note: the section may contain the end/begining of the curve (which is
%closed) meaning that indices alog section may jump to 1
if Peaks(2,1) > Peaks(1,1)
    SecInd(:,1) = Peaks(1,1)+1:Peaks(2,1)-1;
else 
    SecInd(:,1) = [(Peaks(1,1)+1:size(NewDist,1))';(1:Peaks(2,1)-1)'];
end

%If the minimum curvature of the section is greater than 1 - it is
%considered a purely convex section and a common tanget isn't searched. In
%stead the angle of the tangent to the point of min curvatur is taken
if min(K(SecInd)) > 0
    Ang = Theta(SecInd(find(K(SecInd)==min(K(SecInd)))));
    PtInd = SecInd(find(K(SecInd)==min(K(SecInd))));
else
    %If not removes concave (or with curvature < 1) points from the section
    %as these cannot be the ones through which the tangent is passing. The
    %section is then divided into descrete convex sections (Gaps). between
    %each two there is a chance for the commonn tangent to pass.
    
    SecInd = SecInd(find(K(SecInd)>0));
    ComTanInd = inf;
    for i = 1:size(SecInd,1)-round(size(SecInd,1)*0.15)
        for j = i+round(size(SecInd,1)*0.15):size(SecInd,1)
            aPts = (Ypred(SecInd(i))-Ypred(SecInd(j))) ./ (Xpred(SecInd(i))-Xpred(SecInd(j)));
            bPts = Ypred(SecInd(i)) - aPts*Xpred(SecInd(j));
            PtsDist = sqrt((Xpred(SecInd(i))-Xpred(SecInd(j))).^2+((Ypred(SecInd(i))-Ypred(SecInd(j))).^2));
            if (abs(Theta(SecInd(i)) - Theta(SecInd(j))) <= deg2rad(3) ||...
                    abs(Theta(SecInd(i)) - Theta(SecInd(j))) >= deg2rad(177)) &&...
                    abs(mean([Theta(SecInd(i));Theta(SecInd(j))]) - atan(aPts)) <= deg2rad(3)
                if abs(bPts) / PtsDist < ComTanInd
                    ComTanInd = abs(bPts) / PtsDist;
                    Ang = atan(aPts);
                    PtInd = [SecInd(i),SecInd(j)];
                end
            end
        end
    end
end
                    
            
    
    
%                     % if difference between the tangent is less than 1
%                     if abs(AngPt1-AngPt2) <= deg2rad(5) || abs(AngPt1-AngPt2) >= deg2rad(175)
%                         aPt1 = tan(AngPt1);
%                         aPt2 = tan(AngPt2);
%                         bPt1 = Ypred(PGap1) - aPt1*Xpred(PGap1);
%                         bPt2 = Ypred(PGap2) - aPt2*Xpred(PGap2);
%                         aPts = (Ypred(PGap1)-Ypred(PGap2)) ./ (Xpred(PGap1)-Xpred(PGap2));
%                         bPts = Ypred(PGap1) - aPts*Xpred(PGap1);
% %                         if abs(mean([aPt1;aPt2])-aPts) <= tan(deg2rad(3)) %abs(bPt1-bPt2) < 0.01*min(abs([bPt1;bPt2]))
%                             GComTans(m,1:2) = [PGap1, PGap2];
%                             GComTans(m,3) = atan(aPts);
%                             GComTans(m,4) = abs(bPts);
%                             GComTans(m,5) = sqrt((Xpred(PGap1)-Xpred(PGap2)).^2+((Ypred(PGap1)-Ypred(PGap2)).^2));
%                             GComTans(m,6) = abs(mean([aPt1;aPt2])-aPts);%abs(bPt1-bPt2);
%                             %Compare the Y of each point with the
%                             %prediction of its corresponding pt and slope.
%                             %sum and search for the min.a
%                             m=m+1;
% %                         end
%                     end
%                 end
%             end
%            if exist('GComTans','var')
%                [~, Ind] = sort(GComTans(:,6));
%                GComTans = GComTans(Ind,:);
%                if size(GComTans,2)>10
%                    GComTans = GComTans(1:10,:);
%                end
%                ComTans(n,:) = GComTans(find(GComTans(:,4)==min(GComTans(:,4))),:);
%                n=n+1;
%                clear GComTans 
%            end
%         end
%     end
%     ComTans = ComTans(find(ComTans(:,4)==min(ComTans(:,4))),:);
%     Ang = ComTans(1,3);
%     PtInd = ComTans(1,1:2);
% %     ComTans = ComTans(find(ComTans(:,4)>=max(ComTans(:,4))-0.1*max(ComTans(:,4))),:);
% %     Ang = ComTans(find(ComTans(:,5)==max(ComTans(:,5))),3);
% %     PtInd = ComTans(find(ComTans(:,5)==max(ComTans(:,5))),1:2);
% end
            
end

