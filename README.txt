RUNNING THE FINAL MATLAB FILES

1. Open MATLAB and set Current Folder to this folder.
2. Enter this command in the Command Window:

   RUN_PAPER

This computes the two-mode MFC results and exports the final numerical plots.
No previous results, Python installation or other project folders are needed.
Keep all .m files and colormaps.mat together. The MAT file contains only the
colour palettes needed to reproduce the paper's plots.

Required: MATLAB R2026a with Optimization Toolbox (tested with Update 5).
Allow several minutes. The run creates data/, figures/ and verification/
folders for its outputs. Figures are saved as editable .fig, PDF and PNG files.
Rerunning regenerates those outputs. It preserves the existing MATLAB path.

Figure 1 is a separately authored schematic; MATLAB produces Figures 2-10
and C.1-C.4. The model uses the paper's two-mode MFC assumptions and 400 V.
The paper's stated d31 linear-extrapolation limitation still applies.
