% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [Cmax,CDmax,Cmin,CDmin,vN] = SurfPrCurv(Tr)
%SURFPRCURV Returns the principal curvatures and directions based on
%surfasce fitting
%   Based on Curvature Estimation of 3D Point Cloud Surfaces Through the Fitting of Normal
% Section Curvatures ASIAGRAPH 2008 PROCEEDINGS and A Novel Cubic-Order Algorithm for Approximating
% Principal Direction Vectors ACM Transactions on Graphics 2004
%
% Inputs:
%   Tr
% Outputs:
%   Cmax
%   CDmax
%   Cmin
%   CDmin
%   vN
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%Triangulate and calculate face normals, circumcenters, and vertex normals
% tic
f = Tr.ConnectivityList;
v = Tr.Points;
vN = vertexNormal(Tr);

%Neighborhood is found through K nearest neighbors search. Despite using
%Euclidean distance, is high-resolution meshes it provides good results due
%to the small neigborhood represented by many points. Also the least squares
%solution smoothes problems out.
vNeiSz = 30;
Vatt = knnsearch(v,v,'k',vNeiSz);


% assign coords and normals to pages in 3D arrays
% translated and rotated to local coord system of p(1) = 0 0 0, N(1) = 0,0,-1
% zero prelocation for speed
vNei = zeros(min(vNeiSz),3,size(v,1));
NorNei = zeros(min(vNeiSz),3,size(v,1));
Axis = zeros(size(v,1),3);
ang = zeros(size(v,1),1);
for i=1:size(v,1)
    Axis(i,:) = cross(vN(Vatt(i,1),:),[0,0,-1],2);
    ang(i) = acos(dot(vN(Vatt(i,1),:),[0,0,-1],2));
    Rmat = RMatAxAng(Axis(i,:),-ang(i));
    vNei(:,:,i) = (Rmat*(v(Vatt(i,:),:)-v(Vatt(i,1),:))')';
    NorNei(:,:,i) = (Rmat*vN(Vatt(i,:),:)')';
end
%toc


% Least squares problem formulation following https://doi.org/10.1145/966131.966134 
% Note that A and b are in 3D where each page represents a neighborhood and
% is solved using "pagemldevide"
A = [[0.5.*vNei(:,1,:).^2, vNei(:,1,:).*vNei(:,2,:), 0.5.*vNei(:,2,:).^2, vNei(:,1,:).^3,...
        vNei(:,1,:).^2.*vNei(:,2,:), vNei(:,1,:).*vNei(:,2,:).^2, vNei(:,2,:).^3];...
        [vNei(:,1,:), vNei(:,2,:), zeros(size(vNei,1),1,size(vNei,3)), 3.*vNei(:,1,:).^2,...
        2.*vNei(:,1,:).*vNei(:,2,:), vNei(:,2,:).^2, zeros(size(vNei,1),1,size(vNei,3))];...
        [zeros(size(vNei,1),1,size(vNei,3)), vNei(:,1,:), vNei(:,2,:),...
        zeros(size(vNei,1),1,size(vNei,3)), vNei(:,1,:).^2, 2.*vNei(:,1,:).*vNei(:,2,:),...
        3.*vNei(:,2,:).^2]];

b = [vNei(:,3,:); -NorNei(:,1,:)./NorNei(:,3,:); -NorNei(:,2,:)./NorNei(:,3,:)];

Coefs = squeeze(pagemldivide(A,b));
%toc

% eigenvectors and eigenvalues of each Weingarten matrix corrspond to
% magnitude and directions of principal curvatures. Note that directions
% are re-rotated to original coordinate system
% zeros prelocation for speed

Vecs = zeros(2,2,size(v,1));
Vals = zeros(2,2,size(v,1));
Vals2 = zeros(2,size(v,1));
Vecs1 = zeros(3,size(v,1));
Vecs2 = zeros(3,size(v,1));
for i =1:size(v,1)
    [Vecs(:,:,i), Vals(:,:,i)] = eig([Coefs(1,i),Coefs(2,i);Coefs(2,i),Coefs(3,i)]);
    Vals2(1:2,i) = [Vals(1,1,i);Vals(2,2,i)];
    Rmat = RMatAxAng(Axis(i,:),ang(i));
    VecsT = Rmat*[[Vecs(1,:,i),0];[Vecs(2,:,i),0];[0,0,-1]];
    Vecs1(:,i) = VecsT(:,1);
    Vecs2(:,i) = VecsT(:,2);
end
%toc

%remove outliers
% Vals2(1,:) = filloutliers(Vals2(1,:),'clip');
% Vals2(2,:) = filloutliers(Vals2(2,:),'clip');

OutlrsMin = find(isoutlier(Vals2(1,:),"quartiles"));
OutlrsMax = find(isoutlier(Vals2(2,:),"quartiles"));

vM = v;
vM([OutlrsMax,OutlrsMin],:) = inf(size([[OutlrsMax,OutlrsMin],1],3));
VattNoMax = knnsearch(vM,v(OutlrsMax,:),'k',50);
DistMax = zeros(size(VattNoMax));
for i=1:size(OutlrsMax,2)
    DistMax(i,:) = pdist2(v(OutlrsMax(i),:),v(VattNoMax(i,:),:));
end
DistMax = (1./DistMax)./sum((1./DistMax),2);
NewValsMax = reshape(Vals2(2,VattNoMax),size(VattNoMax));
NewValsMax = sum(DistMax.*NewValsMax,2);
Vals2(2,OutlrsMax) = NewValsMax;
Cmax = Vals2(2,:)';
% Cmax(OutlrsMax) = NewValsMax;
CDmax = Vecs2';

VattNoMin = knnsearch(vM,v(OutlrsMin,:),'k',50);

DistMin = zeros(size(VattNoMin));
for i=1:size(OutlrsMin,2)
    DistMin(i,:) = pdist2(v(OutlrsMin(i),:),v(VattNoMin(i,:),:));
end
DistMin = (1./DistMin)./sum((1./DistMin),2);
NewValsMin = reshape(Vals2(2,VattNoMin),size(VattNoMin));
NewValsMin = sum(DistMin.*NewValsMin,2);
Vals2(1,OutlrsMin) = NewValsMin;
Cmin = Vals2(1,:)';
% Cmin(OutlrsMin) = NewValsMin;
CDmin = Vecs1';

% %Split into min and max curvature and flip accordingly
% if max(abs(Vals2(1,:)),[],2) > max(abs(Vals2(2,:)),[],2)
%     Cmax = Vals2(1,:)';
%     CDmax = Vecs1';
%     Cmin = Vals2(2,:)';
%     CDmin = Vecs2';
% else
%     Cmax = Vals2(2,:)';
%     CDmax = Vecs2';
%     Cmin = Vals2(1,:)';
%     CDmin = Vecs1';
% end
% if abs(min(Cmax))>max(Cmax)
%     Cmax = -Cmax;
%     Cmin = -Cmin;
% end
%toc

end

%%

%Old segments using triangulation - too slow for large meshes
% N=neighbors(Tr);
% % toc
% 
% 
% %find vertex attachments and clip to uniform size of minimum cell for array
% %vectorization assignment - Note that placing PiD in first place is
% %conducted here for this reason
% tic
% Vatt = vertexAttachments(Tr);
% % Vatt = cellfun(@(x) unique(f(unique(N(x,:)),:)),Vatt,'UniformOutput',false);
% % for i =1:2
% %     Vatt = cellfun(@cell2mat, cellfun(@(x) vertexAttachments(Tr,x)', Vatt, 'UniformOutput', false), 'UniformOutput', false);
% %     Vatt = cellfun(@(x) unique(f(x,:)), Vatt, 'uniformoutput', false);
% % end
% 
% 
% Vatt = cellfun(@(x) reshape(neighbors(Tr, x'),[],1),Vatt,'uniformoutput',false);
% for i=1:8
%     Vatt = cellfun(@(x) unique(neighbors(Tr, x)),Vatt,'uniformoutput',false);
%     toc
% end
% vNeiSz = cellfun(@(x) numel(x),Vatt);
% Vatt = cellfun(@(x) x(1:min(vNeiSz))', Vatt, 'UniformOutput', false);
% Vatt = cell2mat(Vatt);
% toc
% Vatt = cellfun(@(x) unique(f(x,:)),Vatt,'uniformoutput',false);
% toc
% PiD = num2cell((1:size(v,1))');
% Vatt = cellfun(@(x,y) setdiff(x,y), Vatt, PiD, 'UniformOutput', false);
% Vatt = cellfun(@(x,y) [y;x], Vatt, PiD, 'UniformOutput', false);
% vNeiSz = cellfun(@(x) numel(x),Vatt);
% Vatt = cellfun(@(x) x(1:min(vNeiSz))', Vatt, 'UniformOutput', false);
% Vatt = cell2mat(Vatt);
% toc

% tic
% Vatt = vertexAttachments(Tr);
% Fatt = cellfun(@(x) neighbors(Tr,x'),Vatt,'uniformoutput',false);
% for i=1:5
%     Fatt = cellfun(@(x) neighbors(Tr,unique(x)),Fatt,'uniformoutput',false);
% end
% Vatt = cellfun(@(x) unique(f(reshape(unique(x),[],1),:)),Fatt,'UniformOutput',false);
% vNeiSz = cellfun(@(x) numel(x),Vatt);
% Vatt = cellfun(@(x) x(1:min(vNeiSz)),Vatt,'uniformoutput',false);
% Vatt = [(1:size(Vatt,1))',cellfun(@(x) cell2mat(x))];
% min(vNeiSz)
% toc


