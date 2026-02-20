% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [SegRamV, SegRamF, NchMin] = SegRam(v, f, Side)
%SEGRAM returns the vertices and triangles indexed of either ramus or both
%given the vertices and tringles of a complete positioned mandible and a
%user selection
%   SEGRAM rotates the vertices about y to standard position (on table,
%   symmetry plane ZY, chin towards positive Z). it returns a newly indexed
%   vertices and faces list according to user selection. Note that SegRamV 
%   contains indices for v (which is still not fully rotated in the calling
%   function's workspace) whole f is a triangle list indexing v(SegRamV,:).
%   If both Rami are to be segmented SegRam returns two 2x1 cell arrays
%   with the first contating the right, and the second the left.
%
% Inputs:
%   v
%   f
%   Side
% Outputs:
%   SegRamV
%   SegRamF
%   NchMin
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.


switch Side
    case 'Right'
        [SegRamV, SegRamF,NchMin] = ReInd(v,f,Side);
    case 'Left'
        [SegRamV, SegRamF,NchMin] = ReInd(v,f,Side);
    case 'Both'
        [SegRamRV, SegRamRF,NchMinR] = ReInd(v,f,'Right');
        [SegRamLV, SegRamLF,NchMinL] = ReInd(v,f,'Left');
        SegRamV = {SegRamRV;SegRamLV};
        SegRamF = {SegRamRF;SegRamLF};
        NchMin = {NchMinR;NchMinL};
end
end

function [SegV, SegF, NchMin] = ReInd(v,f,Hem)

%Measures the max "hight" of the positioned right hemimandible at evenly
%spaced bins and differentiates
if strcmp(Hem,'Right')
    Hem = find(v(:,1)<=0);
else
    Hem = find(v(:,1)>=0);
end
ZSmpInd =  linspace(min(v(Hem,3)),max(v(Hem,3)));
ZGrpInd = discretize(v(Hem,3),ZSmpInd);
MndHgt = zeros(1,max(ZGrpInd));
for i =1:max(ZGrpInd)
    MndHgt(i) = max(v(Hem(ZGrpInd==i),2));
end
DzDx = gradient(MndHgt,ZSmpInd(2)-ZSmpInd(1));

%Loop to find the longest "rise" in hemimandible height, representing the
%ramus acesnsion.
SegLen=1;
SegCount=1;
for i =2:size(DzDx,2)
    if DzDx(i)<0
        SegLen = SegLen+1;
    elseif DzDx(i)>0 && DzDx(i-1)<0
        GhtDiff = MndHgt(i-1-SegLen) - MndHgt(i-1);
        LengthPos(SegCount,1:3) = [i,SegLen,GhtDiff];
        SegCount = SegCount+1;
        SegLen=1;
    end
end
Pos = find(LengthPos(:,2).*LengthPos(:,3)==max(LengthPos(:,2).*LengthPos(:,3)));
if numel(Pos)>1
    Pos = Pos(end);
end
NchMin = [MndHgt(LengthPos(Pos-1,1)-1),ZSmpInd(LengthPos(Pos-1,1)-1)];


%Re indexing the v f list
vRamInd = find(v(Hem,3)<=ZSmpInd(LengthPos(Pos,1)-1));
vRamInd = Hem(vRamInd);
fRamInd = find(any(ismember(f,vRamInd),2));
fRam = f(fRamInd,:);
[SegV,~,IndUptN] = unique(fRam);
SegF = reshape(IndUptN,size(fRam));


end

