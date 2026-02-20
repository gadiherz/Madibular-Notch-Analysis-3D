Mandibles_WF_Clean (interactive workflow)

Authors and copyright
- Authored by Gadi Herzlinger and Uzy Smilansky.
- Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
- Documentation assistance: ChatGPT (comments and documentation only).

What this repository contains
- A cleaned, minimal, dependency-trimmed version of the mandibular upper-ramus workflow.
- The computational code is kept intact; cleanup focuses on:
  (i) separating core functions vs. entry scripts, (ii) retaining only required dependencies,
  (iii) adding a simple startup and an explicit entry-point script, and
  (iv) adding comments and lightweight documentation.

Folder structure
- workflow/
    Rak_Line_Extraction_WF.m         Main interactive script: load model, orient, segment, extract.
    RamLineWF.m                      Computes curve-based descriptors from extracted curves.
    run_workflow_interactive.m       Convenience wrapper that runs the two scripts above.
    run_batch_export_excel.m         Comprehensive batch wrapper with an excel output of variables
- src/
    Core functions used by the workflow (curvature estimation, ridge tracing, curve processing).
- apps/
    CompFragSelect.mlapp             App used for interactive completeness/fragment selection.
    SelectSide4Pos.fig/.m            UI used when automatic table-positioning is ambiguous.
- docs/
    EXPECTED_INTERACTION_FLOW.md     Step-by-step description of user prompts/choices.
    PIPELINE_OVERVIEW.md             High-level call graph and data-flow notes.
    DEPENDENCIES.md                  MATLAB/toolbox and runtime assumptions.
    Parameter_Registry.md            Full description of parameters and thresholds 

Quick start
1) Open MATLAB and set the current folder to this repository root.
2) Run:
     startup
   and then:
     run_batch_export_excel or run_workflow_interactive

Inputs
- The workflow expects a triangulated surface mesh (e.g., .ply/.stl/.obj).
- Model loading is performed in Rak_Line_Extraction_WF.m using readSurfaceMesh (availability depends on MATLAB release/toolboxes).

Outputs
- The pipeline is workspace-oriented (consistent with the original implementation).
- Primary per-specimen results are stored in the struct array Res (created/updated by the workflow scripts).
- Downstream analyses (e.g., plotting, export, statistics) are intentionally left outside this repository to preserve original results.

Notes on reproducibility and scope
- The workflow is inherently interactive (Apps/UI decisions are part of the processing). For publication-quality use, record the user choices made during the session (see docs/EXPECTED_INTERACTION_FLOW.md).
- This repository does not change algorithms or parameter values; it only reorganizes files, trims unused dependencies, and adds documentation to make the codebase readable and maintainable.
