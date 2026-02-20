% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [NewDist,Wn,CoefsX,CoefsY,Theta,K,Xpred,Ypred,TotLen] = FourierCoords(v,Neff)
%This script is used to get the coefficient for the X and Y cartesian
%coordinates as function of the cumulative distance along the a closed
%curve, the X and Y of the points used to calculate them. It requries the
%rotated point cloud given by SymmTens_2.m and operates FourierCoefs.m (see
%documentation therein). It also calculates the smoothing weights parameter
%Wn used in the function. It then provides a cumulative distance vector 
%(fixed interval - NewDist) which are positioned in the fuction to predict
%the X and Y's along the curve. This is good to visualize whether the
%function works, but also to differentiate for varios purposes (See
%CurvCalc.m and TanCalc.m). 
%Note that while FourierCoefs.m has been modified to handle 3D points, this
%script is used by the workflow only in 2D cases.
%
% Inputs:
%   v
%   Neff
% Outputs:
%   NewDist
%   Wn
%   CoefsX
%   CoefsY
%   Theta
%   K
%   Xpred
%   Ypred
%   TotLen
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.

clear Xs CoefsX Ys CoefsY TDist CuDist Wn Xpred Ypred NewDist Theta K

%Get the coordinates of points used in calculation their cumulative distance
%along the curve and the coefficients for X and Y (also the total distance
%which is normalized so always equal to 1).

[TotLen, CoefsX, CoefsY, ~, ~] = FourierCoefs(v); %, 0.5);

%Calculate Smoothing weight parameter 1/(1+exp(1+(n-N_effective)/Delta).
%You can play with N_effective and Delta to determine the smoothing of the
%calculated curve with lower N_effective and higher Delta meaning smoother
%but less detailed curve. Neff is given to the function. in cases where no
%smoothing is required it is set to Inf.


Delta = 3;
if Neff==inf
    Wn = ones(size(CoefsX,1),1);
else
    Wn = 1./(1+exp(1+(((1:1:size(CoefsX,1))-Neff)/Delta)))';
end


% create vector of cumulative distances 800 points 0->1
NewDist = (0:TotLen/800:TotLen)';

%position in formula to get X and Y (predicted). Note that the X and Y
%center (of the original points used to calculate the coefficients) needs
%to be added to the sum of terms provided by Sum4Fourier (see documentation
%therein). Again, this script is only used for the 2D cases.
for i = 1:size(NewDist,1)
    Xpred(i,1) = Sum4Fourier(TotLen,CoefsX,NewDist(i),Wn);
    Ypred(i,1) = Sum4Fourier(TotLen,CoefsY,NewDist(i),Wn);
end

%Calculate the angle between the tangent to each point on the curve (theta)
%and the X axis. This will allow the subsequent rotation of the point cloud
%(and later the whole model). See documentation in TanCalc.m
%This is in fact the first derivative so that theta(s)=arctan(y'(s)/x'(s)).

Theta = TanCalc(NewDist,Wn,CoefsX,CoefsY);

%Because tan is deifned over -pi <-> pi the result has to be manipulated
%by adding pi to one half and 2*pi to the second hlaf (the for loop of
%Th1). These variuos types of manipulation are not used anymore because
%they just complicate matters and can be worked without. Code saved as 
%comment for potential necessities.

%Theta = Theta - (0.5*pi);

% Th1 = Theta(1);
% m=0;
% for i=2:size(Theta,1)
%     if Theta(i)-Theta(i-1)>2.95
%         m=m-1;
%     elseif Theta(i)-Theta(i-1)<-2.95
%         m=m+1;
%     end
%     Th1(i,1) = Theta(i)+(pi*m);
% end
% if abs(Theta(end)-Theta(1))>2.95
%     Th1(1,1) = Theta(1)+(pi*(m+1));
% end
% Theta = Th1;
% clear Th1

%Calculate the curvature at each point along the curve. See documentation 
%in CurvCalc.m This is in fact the second derivative so that 
%K(s) = Theta'(s) = (y'(s)*x''(s) -x'(s)*y''(s))/x'(s)+y'(s)
K = CurvCalc(NewDist,Wn,CoefsX,CoefsY);


%Plot the original curve (Xs, Ys), and the predicted points and curve
%(Xpred, Ypred) and the theta and curvature as functions of the cumulative 
%distance at fixed intervals 0->1 (NewDist).
%XYKTPlotter; %- Stop for now (incorporation into main workflow in
%SymmTens_2.m)
end