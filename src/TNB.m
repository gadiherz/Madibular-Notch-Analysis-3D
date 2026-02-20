% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function [T,N,B,kappa,tau] = TNB(x,y,z,h)
%
% Inputs:
%   x
%   y
%   z
%   h
% Outputs:
%   T
%   N
%   B
%   kappa
%   tau
%
% Notes:
%   - Documentation-only additions; computations unchanged.
%   - See the Supplementary Online Methods (SOM) for mathematical details.
	
	if(nargin == 2)
		z = zeros(size(x));
	end
	
	dr = [gradient(x(:),h) , gradient(y(:),h) , gradient(z(:),h)];
	ds = sum(dr.^2,2).^(0.5); % Arc length associated with each point. ||dr||.
	
	T = dr ./ ds; % Unit tangent vector.
	
	dT = [gradient(T(:,1),h) , gradient(T(:,2),h) , gradient(T(:,3),h)]; % T'(t).
	dTds = dT ./ ds;
	
	kappa = sum((dT./ds).^2,2).^(0.5);
	
	N = dTds ./ kappa; % Unit normal vector.
	
	B = cross(T,N); % Unit bi-normal vector.
	
	dB = [gradient(B(:,1),h) , gradient(B(:,2),h) , gradient(B(:,3),h)];
	
	tau = dot(-dB./ds,N,2); % Torsion.
end