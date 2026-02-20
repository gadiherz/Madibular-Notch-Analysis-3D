EXPECTED_INTERACTION_FLOW

This workflow is interactive. The steps below describe the prompts and user actions that may appear during a typical run.

Entry points
- Preferred: run_workflow_interactive.m
- Equivalent manual sequence:
  1) startup.m
  2) Rak_Line_Extraction_WF.m
  3) RamLineWF.m

Stage 1: Rak_Line_Extraction_WF.m (interactive extraction)
1) File selection dialog
   - A file dialog prompts you to select a mesh file (*.ply / *.stl / *.obj).

2) CompFragSelect (App)
   The App opens and asks:
   - “Is the input model a complete mandible or ramus fragment?”
     Options: Complete / Fragment

   A) If “Complete”
      - The workflow attempts automated positioning (“PosOnTable”).
      - A prompt asks: “Which sides can be used?”
        Options: Both / Left / Right
      - Ramus segmentation and initial landmarking run automatically.

   B) If “Fragment”
      - A prompt asks: “Is the input model a left or right fragment?”
        Options: Left / Right
      - In the UIAxes, candidate curvature peaks are shown.
      - You click twice:
        (i) “Select the middle of the mandibular notch”
        (ii) “Select the condyle”
      - A confirmation dialog appears:
        “Is the positioning OK?”
        Options: OK / Restart
        Choosing Restart returns to the initial selection stage.

3) Optional: SelectSide4Pos (GUI)
   - In some cases, an additional GUI may appear if the automatic positioning/orientation is ambiguous.
   - Follow the on-screen selection to choose the correct side/orientation.

Stage 2: RamLineWF.m (descriptor computation)
- No user interaction is expected.
- The script assumes the workspace variables created by Rak_Line_Extraction_WF.m are available in the same session.

Outputs and what to save
- The workflow leaves results in the MATLAB workspace.
- The primary container is the struct array Res (per specimen/side), which is updated across stages.
- For reproducibility, save:
  - the workspace (or at minimum Res and the extracted curve vertex indices),
  - the selected file name(s),
  - the interactive decisions (Complete/Fragment; side; Restart/OK choices).
