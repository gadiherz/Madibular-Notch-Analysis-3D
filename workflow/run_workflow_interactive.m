% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------

% run_workflow_interactive.m
% Entry point for the interactive mandibular upper-ramus workflow.
%
% This wrapper:
%   - ensures paths are set (calls startup.m)
%   - runs the two legacy scripts in sequence:
%       1) Rak_Line_Extraction_WF.m (interactive positioning/selection + feature extraction)
%       2) RamLineWF.m              (curve resampling + differential-geometry descriptors)
%
% Outputs are left in the workspace, consistent with the original pipeline.
%
% Important:
%   The underlying workflow uses interactive steps (e.g., App-based selection).
%   This script is intended as a convenience wrapper and does not alter results.

% Dependencies:
%   - startup.m (adds required paths)
%   - apps/CompFragSelect.mlapp and apps/SelectSide4Pos.fig/.m
%   - src/ functions (added to path by startup)
%   - readSurfaceMesh (mesh import; availability depends on MATLAB release/toolboxes)
%
startup;

% Step 1: load + orient + segment + extract candidate curves/surfaces
Rak_Line_Extraction_WF;

% If the user cancelled selection / extraction, Rak_Line_Extraction_WF returns
% without producing Res. In that case, stop cleanly.
if ~exist('Res','var') || isempty(Res)
    return
end

% Step 2: resample curves + compute TNB frames, curvature/torsion etc.
RamLineWF;
