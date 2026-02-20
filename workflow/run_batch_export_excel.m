% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% run_batch_export_excel.m
%
% Purpose:
%   Run the existing workflow across multiple mesh files, collect the 12
%   variables used in PCA/JMP (one row per analyzed ramus), and export them
%   to an Excel file.
%
% Notes:
%   - This wrapper does not modify extraction computations.
%   - The per-specimen workflow remains interactive (CompFragSelect).
%   - All accumulator variables are prefixed with MandiblesWF_ so they
%     survive the workflow's batch-safe clearvars step.
% -------------------------------------------------------------------------

startup;

% Determinism / repeatability
rng(0,'twister');
try
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
catch
    % If Parallel Computing Toolbox is unavailable, ignore.
end

MandiblesWF_BatchMode  = true;
MandiblesWF_NoFigures  = true;

% Select input meshes
[MandiblesWF_Files, MandiblesWF_Path] = uigetfile({...
    '*.ply;*.stl;*.obj','3D mesh files (*.ply, *.stl, *.obj)';...
    '*.*','All files'}, ...
    'Select 3D model(s) for batch export', ...
    'MultiSelect', 'on');

if isequal(MandiblesWF_Files,0)
    return
end

if ischar(MandiblesWF_Files)
    MandiblesWF_Files = {MandiblesWF_Files};
end

MandiblesWF_FileList = cellfun(@(x) fullfile(MandiblesWF_Path,x), MandiblesWF_Files, 'UniformOutput', false);
MandiblesWF_nFiles   = numel(MandiblesWF_FileList);

% Preallocate results table (metadata + 12 variables)
MandiblesWF_Results = table();

MandiblesWF_i = 1;
while MandiblesWF_i <= MandiblesWF_nFiles

    MandiblesWF_CurrentFile = MandiblesWF_FileList{MandiblesWF_i};

    % Specimen ID from file name
    [~, MandiblesWF_SpecimenID, ~] = fileparts(MandiblesWF_CurrentFile);
    MandiblesWF_SpecimenID = string(MandiblesWF_SpecimenID);

    try
        % Step 1: extraction (interactive)
        Rak_Line_Extraction_WF;

        if ~exist('Res','var') || isempty(Res)
            % No output (e.g., user cancelled)
            MandiblesWF_i = MandiblesWF_i + 1;
            continue
        end

        % Step 2: descriptors
        RamLineWF;

        % --- Extra variables (computed once per mandible) ---
        MandiblesWF_ExpectVolume = exist('Side','var') && iscell(Side);

        MandiblesWF_MandibleVolume = NaN;
        MandiblesWF_NotchLen_R = NaN; MandiblesWF_NotchLen_L = NaN;
        MandiblesWF_CondLen_R  = NaN; MandiblesWF_CondLen_L  = NaN;

        if exist('NRTotLen','var'), MandiblesWF_NotchLen_R = double(NRTotLen); end
        if exist('NLTotLen','var'), MandiblesWF_NotchLen_L = double(NLTotLen); end
        if exist('CRTotLen','var'), MandiblesWF_CondLen_R  = double(CRTotLen); end
        if exist('CLTotLen','var'), MandiblesWF_CondLen_L  = double(CLTotLen); end

        % Mandible volume only for COMPLETE mandibles (Side is cell for complete)
        if MandiblesWF_ExpectVolume && exist('v','var') && exist('f','var') && ~isempty(v) && ~isempty(f)
            MandiblesWF_f = f;
            if min(MandiblesWF_f(:)) == 0
                MandiblesWF_f = MandiblesWF_f + 1; % safety for 0-based faces
            end
            MandiblesWF_p1 = v(MandiblesWF_f(:,1),:);
            MandiblesWF_p2 = v(MandiblesWF_f(:,2),:);
            MandiblesWF_p3 = v(MandiblesWF_f(:,3),:);
            MandiblesWF_MandibleVolume = abs(sum(dot(MandiblesWF_p1, cross(MandiblesWF_p2, MandiblesWF_p3, 2), 2)) / 6);
        end


        % Collect rows (one per side, excluding synthetic N/A row)
        T = struct2table(Res);
        if ismember('Side', T.Properties.VariableNames)
            MandiblesWF_valid = ~(string(T.Side) == "N/A");
        else
            MandiblesWF_valid = true(height(T),1);
        end

        MandiblesWF_r = 1;
        while MandiblesWF_r <= height(T)
            if ~MandiblesWF_valid(MandiblesWF_r)
                MandiblesWF_r = MandiblesWF_r + 1;
                continue
            end

            % Extract variables with robust missing handling
            MandiblesWF_Flags = strings(0,1);

            % Helper: get scalar numeric field or NaN + flag
            MandiblesWF_get = @(name) local_get_scalar(T, MandiblesWF_r, name);

            [TriCrvMean, f1]        = MandiblesWF_get('TriCrvMean');
            [CondyleSurfaceArea, f13] = MandiblesWF_get('CondSurf');
            [CondCurvInt, f2]       = MandiblesWF_get('CondCurvInt');
            [CondPlnAng, f3]        = MandiblesWF_get('CondPlnAng');
            [CondTorInt, f4]        = MandiblesWF_get('CondTorInt');
            [NotchTorInt, f5]       = MandiblesWF_get('NotchTorInt');
            [CondAbsTorInt, f6]     = MandiblesWF_get('CondAbsTorInt');
            [NotchCurvInt, f7]      = MandiblesWF_get('NotchCurvInt');
            [NotchAbsTorInt, f8]    = MandiblesWF_get('NotchAbsTorInt');
            [CondPlnRMSE, f9]       = MandiblesWF_get('CondPlnRMSE');
            [CondRamSymAng, f10]    = MandiblesWF_get('CondRamSymAng');
            [NotchRamSymDev, f11]   = MandiblesWF_get('NotchRamSymDev');
            [RamMandSymAng, f12]    = MandiblesWF_get('RamMandSymAng');

            MandiblesWF_Flags = [MandiblesWF_Flags; f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12; f13];
            MandiblesWF_Flags = MandiblesWF_Flags(MandiblesWF_Flags ~= "");

            if ismember('Side', T.Properties.VariableNames)
                MandiblesWF_Side = string(T.Side(MandiblesWF_r));
            else
                MandiblesWF_Side = "";
            end

                        % --- Extra variables per row (select by side) ---
            NotchCurveLength = NaN;
            CondylarDelineationCurveLength = NaN;

            if MandiblesWF_Side == "Right"
                NotchCurveLength = MandiblesWF_NotchLen_R;
                CondylarDelineationCurveLength = MandiblesWF_CondLen_R;
            elseif MandiblesWF_Side == "Left"
                NotchCurveLength = MandiblesWF_NotchLen_L;
                CondylarDelineationCurveLength = MandiblesWF_CondLen_L;
            else
                MandiblesWF_Flags = [MandiblesWF_Flags; "UNKNOWN_SIDE_FOR_LENGTHS"];
            end

            if isnan(NotchCurveLength)
                MandiblesWF_Flags = [MandiblesWF_Flags; "MISSING_NotchCurveLength"];
            end
            if isnan(CondylarDelineationCurveLength)
                MandiblesWF_Flags = [MandiblesWF_Flags; "MISSING_CondCurveLength"];
            end

            MandibleVolume = MandiblesWF_MandibleVolume;
            if MandiblesWF_ExpectVolume && isnan(MandibleVolume)
                MandiblesWF_Flags = [MandiblesWF_Flags; "MISSING_MandibleVolume"];
            end


            if isempty(MandiblesWF_Flags)
                MandiblesWF_FlagsStr = "";
            else
                MandiblesWF_FlagsStr = strjoin(MandiblesWF_Flags, "; ");
            end

            MandiblesWF_NewRow = table(...
                MandiblesWF_SpecimenID, MandiblesWF_Side, MandiblesWF_FlagsStr, ...
                TriCrvMean, CondCurvInt, CondPlnAng, CondTorInt, NotchTorInt, ...
                CondAbsTorInt, NotchCurvInt, NotchAbsTorInt, CondPlnRMSE, ...
                CondRamSymAng, NotchRamSymDev, RamMandSymAng, MandibleVolume, ...
                CondyleSurfaceArea, NotchCurveLength, CondylarDelineationCurveLength, ...
                'VariableNames', { ...
                    'SpecimenID','Side','Flags', ...
                    'CondyleMeanCurvature_Surface', ...
                    'CondylarDelineationCurvature', ...
                    'CondylarPlaneAngle', ...
                    'CondylarDelineationTorsion', ...
                    'NotchTorsion', ...
                    'CondylarDelineationAbsTorsion', ...
                    'NotchCurvature', ...
                    'NotchAbsTorsion', ...
                    'CondylarDelineatingCurveFitError', ...
                    'CondyleRamusAngle', ...
                    'NotchCurveFitError', ...
                    'MandibleRamusAngle' ...
                    'MandibleVolume', ...
                    'CondyleSurfaceArea', ...
                    'NotchCurveLength', ...
                    'CondylarDelineationCurveLength' ...

                });

            MandiblesWF_Results = [MandiblesWF_Results; MandiblesWF_NewRow];

            MandiblesWF_r = MandiblesWF_r + 1;
        end

    catch ME
        % If anything fails, record a single row with an error flag
        MandiblesWF_Results = [MandiblesWF_Results; table(...
            MandiblesWF_SpecimenID, "", "ERROR: " + string(ME.message), ...
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
            NaN, NaN, NaN, NaN, ...
            'VariableNames', { ...
                'SpecimenID','Side','Flags', ...
                'CondyleMeanCurvature_Surface', ...
                'CondylarDelineationCurvature', ...
                'CondylarPlaneAngle', ...
                'CondylarDelineationTorsion', ...
                'NotchTorsion', ...
                'CondylarDelineationAbsTorsion', ...
                'NotchCurvature', ...
                'NotchAbsTorsion', ...
                'CondylarDelineatingCurveFitError', ...
                'CondyleRamusAngle', ...
                'NotchCurveFitError', ...
                'MandibleRamusAngle' ...
                'MandibleVolume', ...
                'CondyleSurfaceArea', ...
                'NotchCurveLength', ...
                'CondylarDelineationCurveLength' ...
            } )];
    end

    MandiblesWF_i = MandiblesWF_i + 1;
end

% Save Excel
[MandiblesWF_OutFile, MandiblesWF_OutPath] = uiputfile('Results.xlsx','Save results table as', fullfile(MandiblesWF_Path,'Results.xlsx'));
if isequal(MandiblesWF_OutFile,0)
    return
end

MandiblesWF_OutFull = fullfile(MandiblesWF_OutPath, MandiblesWF_OutFile);

try
    writetable(MandiblesWF_Results, MandiblesWF_OutFull);
catch
    % Fallback: write as CSV if Excel write fails
    MandiblesWF_OutFull = regexprep(MandiblesWF_OutFull,'\.xlsx$','.csv','ignorecase');
    writetable(MandiblesWF_Results, MandiblesWF_OutFull);
end

% Also save a MAT copy alongside the spreadsheet
try
    MandiblesWF_MatFull = regexprep(MandiblesWF_OutFull,'\.(xlsx|csv)$','.mat','ignorecase');
    MandiblesWF_T = MandiblesWF_Results; %#ok<NASGU>
    save(MandiblesWF_MatFull, 'MandiblesWF_T');
catch
end


% -------------------------------------------------------------------------
% Local helper: fetch scalar numeric field from a struct2table row
% -------------------------------------------------------------------------
function [val, flag] = local_get_scalar(T, rowIdx, fieldName)
val = NaN;
flag = "";

if ~ismember(fieldName, T.Properties.VariableNames)
    flag = "MISSING_" + string(fieldName);
    return
end

x = T.(fieldName)(rowIdx);

% Allow numeric scalars and 1x1 arrays; everything else becomes NaN
try
    if isempty(x)
        flag = "EMPTY_" + string(fieldName);
        return
    end

    if iscell(x)
        x = x{1};
    end

    if isnumeric(x) && isscalar(x)
        val = double(x);
    else
        % Attempt to convert string-like numerics
        if (isstring(x) || ischar(x))
            tmp = str2double(x);
            if ~isnan(tmp)
                val = tmp;
            else
                flag = "NAN_" + string(fieldName);
            end
        else
            flag = "NONSCALAR_" + string(fieldName);
        end
    end

    if isnan(val)
        flag = "NAN_" + string(fieldName);
    end
catch
    flag = "ERROR_" + string(fieldName);
end
end
