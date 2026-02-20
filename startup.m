% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------

% startup.m
% Adds Mandibles workflow folders to the MATLAB path.
%
% Usage:
%   1) cd into the repository root
%   2) run startup
%
% Notes:
%   - This repository is designed for interactive use (Apps/UI selection).
%   - No computational changes were made relative to the archived workflow;
%     this file only sets paths.

repoRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repoRoot,'src')));
addpath(genpath(fullfile(repoRoot,'workflow')));
addpath(genpath(fullfile(repoRoot,'apps')));
