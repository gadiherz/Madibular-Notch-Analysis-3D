DEPENDENCIES

MATLAB runtime
- This repository is MATLAB-centric and includes:
  - .m functions and scripts
  - a MATLAB App Designer app (.mlapp)
  - a GUIDE GUI (.fig/.m)

Key runtime expectations
- A MATLAB installation that provides:
  - patch/triangulation/graphics primitives
  - readSurfaceMesh for mesh import (availability depends on MATLAB release/toolboxes)
  - basic linear algebra and geometry functions used throughout the src/ code

Optional/toolbox-sensitive functions
- Some MATLAB functions used in this workflow may require specific toolboxes depending on MATLAB version
  (for example, functions related to mesh import, smoothing, or specialized geometry utilities).
- If a function is missing at runtime, MATLAB will report it; resolve by:
  - using a MATLAB release that includes the function, or
  - installing the required toolbox, or
  - swapping import to an equivalent reader (note: algorithm changes are outside the scope of this cleaned package).

External dependencies
- None required for the cleaned package.
- The workflow intentionally avoids hard-coded paths and is designed to run after calling startup.m.

Reproducibility note
- Mesh pre-processing choices (decimation, smoothing) can affect curvature-derived measures.
- In this repository, these steps are preserved as implemented in the legacy workflow; documentation is provided, but no numerical behavior is altered.
