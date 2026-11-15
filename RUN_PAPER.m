function RUN_PAPER()
%RUN_PAPER Compute the current paper results and export final numerical figures.
%   Run from this folder. No saved results or external project files are needed.
root=fileparts(mfilename('fullpath'));
oldPath=path;
cleanup=onCleanup(@()path(oldPath));
addpath(root,'-begin');
assert(strcmp(which('mfc_material'),fullfile(root,'mfc_material.m')));
run_mfc_port(false);
mfc_plot_all(root);
revise_review_figures();
revise_legends();
enlarge_fonts();
fprintf('Finished. Results: %s\nFigures: %s\n',fullfile(root,'data'),fullfile(root,'figures'));
end
