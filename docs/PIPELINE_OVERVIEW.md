PIPELINE_OVERVIEW

Goal
Extract anatomically defined curves/surfaces from a 3D mandibular mesh and compute intrinsic geometric descriptors (curvature/torsion) and plane-relation measures as quantitative variables.

High-level stages and primary code units

0) Startup and path setup
- startup.m
  Adds src/ and apps/ to the MATLAB path.

1) Load model and reduce mesh
- Rak_Line_Extraction_WF.m
  - loads a mesh via readSurfaceMesh
  - performs decimation steps used for curvature calculations

2) Orientation and ramus segmentation
- CompFragSelect (apps/CompFragSelect.mlapp)
  - Complete mandible path:
    - PosOnTable
    - SegRam
  - Fragment path:
    - interactive notch/condyle picking
    - automated alignment to a quasi-anatomical pose

3) Curvature estimation on triangulated surfaces
- SurfPrCurv (src/)
  - computes principal curvature quantities used for ridge/feature extraction

4) Anatomical point detection and condyle extraction
- CrndNchMin, CrndCoron, CondyleSeg2, CondyleSeg3 (src/)
  - detect coronion and notch minimum using geometric criteria
  - isolate condylar region by geometric thresholds

5) Curve extraction
- NotchLine2/NotchLine3 and related helpers (src/)
  - extract:
    - upper condylar surface delineation curve
    - mandibular notch ridge curve

6) Differential-geometry descriptors on extracted curves
- RamLineWF.m (workflow/)
  - curve resampling/smoothing
  - TNB frames
  - curvature and torsion summaries
  - plane association and projection error measures

7) Result packaging
- Res struct array is updated in-place and remains in the workspace.

Design constraints
- The workflow is intentionally script/workspace-driven to preserve legacy results.
- Interactivity is part of the method; record user choices when running publication analyses.
