% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [TotLen, CoefsX, CoefsY, CoefsZ, GCoefs] = FourierCoefs(Bcoords, OpenCurv)
%FourierCoefs finds the coefficients for building the
%Fourier series of the 2D or 3D coordinates as function of arc length.
%
% Inputs:
%   Bcoords
%   OpenCurv
% Outputs:
%   TotLen
%   CoefsX
%   CoefsY
%   CoefsZ
%   GCoefs
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

%It receives an *ordered*(!) list of 2D or 3D coordiates along a closed 
%line. It returns the total length of the closed curve, the XY and Z (when
%applicable) coordinates of the points after sampling and sliding so that
%they are in uniform distance from each other and the coefficients for
%calculating the XY and Z (when applicable) coordinates as a function of
%distance along the curve (normalized 0 to 1) using the Sum4Fourier
%function. The coefficitns for each simension are given in the form of nX2
%so that the first column is the for the Sine and the second for the
%Cosine.

% If not specified curve is assumed to be closed. Note that this is a user
% or code definition as the next section closes the curve formaly (i.e
% Bcoords(:,end) == Bcoords(:,1))
if ~exist('OpenCurv','var')
    OpenCurv = false;
end

%Gets equidistantly sampled points along the curve represented by Bcoords
%and their cumulative distances along it.
[B2, CuDist] = EquiDistPT(Bcoords, OpenCurv);


%This loop is only applied to closed 2D curves and it is used to set a
%constant start point with respect to the coordinate system by sliding.

if size(Bcoords,2) == 2 && ~OpenCurv
    %Centers curve
    B3 = B2 - mean(B2);
    %finds the position on the curve where x changes sign
    PotPts = [];
    if (B3(end,1) < 0 && B3(1,1) > 0) || (B3(end,1) > 0 && B3(1,1)<0) 
        PotPts = [size(B3,1),1];
    end
    iPotPts = [find(any([all([B3(2:end,1) < 0, B3(1:end-1,1) > 0],2),all([B3(2:end,1) > 0, B3(1:end-1,1) < 0],2)],2)),...
        find(any([all([B3(2:end,1) < 0, B3(1:end-1,1) > 0],2),all([B3(2:end,1) > 0, B3(1:end-1,1) < 0],2)],2))+1];
    PotPts = [PotPts; iPotPts];
    RC = find(mean(reshape(B3(PotPts,2),size(PotPts)),2)==min(mean(reshape(B3(PotPts,2),size(PotPts)),2)));
    %Shifts and flips the curve as needed so that the first point will
    %always be the one after x changed sign clocwize
    if B3(PotPts(RC,1),1) < 0 && B3(PotPts(RC,2),1) > 0
        B3 = circshift(B3,-PotPts(RC,2));
    else
        B3 = flip(circshift(B3,-PotPts(RC,2)));
    end
    %returns the curve to place
    B2 = B3 + mean(B2);
    clear  PotPts iPotPts RC B3
end

if OpenCurv
    %The openess of the curve means that it is not periodically smooth
    %which has implication on higher derivatives of the Fourier series.
    %Therefore in open cases a polynomial function is defined so that it is
    %reduced from the orginal curve and the difference between the two
    %function is a periodically smooth function. 
    %Note that it is assomed that open curves are always 3D
    %Mathematical documentation for this part is available seperately
    
    %Defines start and end values for the polynomial and its first and second
    %derivaties. Note that from here the original curve is trimmed by one
    %sampling point in the begining and the end because of the need to
    %calculate the derivative values at the begginnig and end numerically
    CuDist = CuDist(1:end-2);

    
    F0 = [B2(2,1), B2(2,2), B2(2,3)];
    FL = [B2(end-1,1), B2(end-1,2), B2(end-1,3)];
    Ft0 = [(B2(3,1)-B2(1,1))./(2.*CuDist(2)),(B2(3,2)-B2(1,2))./(2.*CuDist(2)),(B2(3,3)-B2(1,3))./(2.*CuDist(2))];
    FtL = [(B2(end,1)-B2(end-2,1))./(2.*CuDist(2)),(B2(end,2)-B2(end-2,2))./(2.*CuDist(2)),...
        (B2(end,3)-B2(end-2,3))./(2.*CuDist(2))];
    Ftt0 = [(B2(3,1)-2.*B2(2,1)+B2(1,1))./CuDist(2).^2,(B2(3,2)-2.*B2(2,2)+B2(1,2))./CuDist(2).^2,(B2(3,3)-2.*B2(2,3)+B2(1,3))./CuDist(2).^2];
    FttL = [(B2(end-2,1)-2.*B2(end-1,1)+B2(end,1))./CuDist(2).^2,...
        (B2(end-2,2)-2.*B2(end-1,2)+B2(end,2))./CuDist(2).^2,...
        (B2(end-2,3)-2.*B2(end-1,3)+B2(end,3))./CuDist(2).^2,];
    
    %uses a matrix to solve a system of linear equations to get the
    %coeffitients of the polynomial
    
    AMat = [1,CuDist(1),CuDist(1).^2,CuDist(1).^3,CuDist(1).^4,CuDist(1).^5;...
        0,1,2*CuDist(1),3*CuDist(1).^2,4*CuDist(1).^3,5.*CuDist(1).^4;...
        0,0,2,6.*CuDist(1),12.^CuDist(1).^2,20.*CuDist(1).^3;...
        1,CuDist(end),CuDist(end).^2,CuDist(end).^3,CuDist(end).^4,CuDist(end).^5;...
        0,1,2*CuDist(end),3*CuDist(end).^2,4*CuDist(end).^3,5.*CuDist(end).^4;...
        0,0,2,6.*CuDist(end),12.*CuDist(end).^2,20.*CuDist(end).^3];
    
    RVec = [F0;Ft0;Ftt0;FL;FtL;FttL];
    GCoX = AMat\RVec(:,1);
    GCoY = AMat\RVec(:,2);
    GCoZ = AMat\RVec(:,3);
    
    %GCoefs contains six coeffs for each dimension of the curve
    GCoefs = [GCoX,GCoY,GCoZ];
    
    %samples points for the polynomial in accordance with the original
    %sampling of the curve (CuDist)
    GXs = GCoX(1) + GCoX(2)*CuDist + GCoX(3)*CuDist.^2 + GCoX(4)*CuDist.^3 + GCoX(5)*CuDist.^4 + GCoX(6)*CuDist.^5;
    GYs = GCoY(1) + GCoY(2)*CuDist + GCoY(3)*CuDist.^2 + GCoY(4)*CuDist.^3 + GCoY(5)*CuDist.^4 + GCoY(6)*CuDist.^5;
    GZs = GCoZ(1) + GCoZ(2)*CuDist + GCoZ(3)*CuDist.^2 + GCoZ(4)*CuDist.^3 + GCoZ(5)*CuDist.^4 + GCoZ(6)*CuDist.^5;
    clear F0 FL Ft0 FtL Ftt0 FttL AMat RVec GCoX GCoY GCoZ

    %     %For testing new function
%     Fx = B2(2:end-1,:) - [GXs,GYs,GZs];
%     
% %     For testing first derivative
%     GtXs = GCoefs(2,1)+ 2.*GCoefs(3,1)*CuDist + 3.*GCoefs(4,1)*CuDist.^2 + 4.*GCoefs(5,1)*CuDist.^3 + 5.*GCoefs(6,1)*CuDist.^4;
%     GtYs = GCoefs(2,2)+ 2.*GCoefs(3,2)*CuDist + 3.*GCoefs(4,2)*CuDist.^2 + 4.*GCoefs(5,2)*CuDist.^3 + 5.*GCoefs(6,2)*CuDist.^4;
%     GtZs = GCoefs(2,3)+ 2.*GCoefs(3,3)*CuDist + 3.*GCoefs(4,3)*CuDist.^2 + 4.*GCoefs(5,3)*CuDist.^3 + 5.*GCoefs(6,3)*CuDist.^4;
%     
%     B2tN1 = (B2(3:end,:) - B2(1:end-2,:))./(2.*CuDist(2));
%     B2tN1 = B2tN1./vecnorm(B2tN1,2,2);
%     Ftx = B2tN1 - [GtXs,GtYs,GtZs];
    
    %CuDist used to calculate Gs has already been trimmed for the
    %polynomial and now it has to be trimmed again but only at the end for
    %the Fourier series
    Xs = B2(2:end-2,1)-GXs(1:end-1);
    Ys = B2(2:end-2,2)-GYs(1:end-1);
    if size(Bcoords,2) == 3
        Zs = B2(2:end-2,3)-GZs(1:end-1);
    end
    %clear GXs GYs GZs
else
    %If curve is closed it can simply be used after the removal of the
    %duplicacy by removing the last point
    Xs = B2(1:end-1,1);
    Ys = B2(1:end-1,2);
    if size(Bcoords,2) == 3
        Zs = B2(1:end-1,3);
    end
    %Return NaNs if curve is closed
    GCoefs = NaN(6,3);
end

%The total length of the curve (after adustments and trimming) is kept for
%returning and CuDist is trimmed by one sample point at the end to avoid
%duplicacy (regardless of openess or dimensionality).
TotLen = CuDist(end);
CuDist = CuDist(1:end-1);


%The Fourier coefficients are found by solving a system of linear
%equations, namely fitting the points on the curve (original if closed,
%modified if open). CuDist runs from 0 to almost the end of the curve (without
%duplicacy). The number of coefficients equals the number of points used
%diveded by two (must be even). The first coefficients (namely n=0 where the
%sine term equals 0) are included in the coefficients and not calculated
%seperately. Note that the matrix is inverted using the pseudoinverse
%function to avoid singularity cause by inclusion of the first terms. The
%odd colums in the matrix are for the cosine and the even ones for the
%sines

Mat = zeros(length(Xs),size(Xs,1));
for i=1:(size(Xs,1)-1)/2
    Mat(:,2*i) = sin(2*pi*i*CuDist/TotLen);
    Mat(:,2*i+1) = cos(2*pi*i*CuDist/TotLen);
end
Mat(:,1) = ones(length(Xs),1);
CoefsX = inv(Mat)*Xs;
CoefsY = inv(Mat)*Ys;
if size(Bcoords,2) == 3
    CoefsZ = inv(Mat)*Zs;
    CoefsZ = [[0;CoefsZ(2:2:end)],CoefsZ(1:2:end)];
end

CoefsX = [[0;CoefsX(2:2:end)],CoefsX(1:2:end)];
CoefsY = [[0;CoefsY(2:2:end)],CoefsY(1:2:end)];
%Return NaNs for Z if curve is 2D

if size(Bcoords,2) ~= 3
    Zs = NaN;
    CoefsZ = ones(size(CoefsX)) * NaN;
end


end


