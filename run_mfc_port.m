function results = run_mfc_port(makePlots)
%RUN_MFC_PORT Compute the current MFC model entirely in MATLAB.
%   RUN_MFC_PORT(false) computes all numerical results.
%   RUN_MFC_PORT(true) also renders the publication figure set.
%   The calculation never reads saved reference results.
%   Use RUN_PAPER to reproduce the final manuscript plot presentation.
if nargin==0,makePlots=true;end
here=fileparts(mfilename('fullpath'));root=here;
oldPath=path;restore=onCleanup(@()path(oldPath));addpath(here,'-begin');
data=fullfile(root,'data');if ~isfolder(data),mkdir(data);end
checksDir=fullfile(root,'verification');if ~isfolder(checksDir),mkdir(checksDir);end
started=tic;
[models,material]=mfc_material;
mfc_write_json(fullfile(data,'source_material.json'),material);
mfc_write_json(fullfile(data,'calculated_models.json'),models);
[results,summary,benchmarks,losses,convergence]=mfc_analyse(models);
[results,extensionSummary,extensionBenchmarks,angleRows,frontRows,modeComparison]=mfc_extensions(results);
names={'d33','d31'};
for k=1:2
    name=names{k};atlas=results.(name).atlas;extensions=results.(name).extensions;
    save(fullfile(data,[name '_atlas.mat']),'-struct','atlas','-v7');
    save(fullfile(data,[name '_design_extensions.mat']),'-struct','extensions','-v7');
    tableData=table(flat(atlas.phi),flat(atlas.psi),flat(atlas.Ex)/1e9,flat(atlas.Ey)/1e9, ...
        flat(atlas.Gel)/1e9,flat(atlas.G(:,:,1)),flat(atlas.G(:,:,2)), ...
        flat(atlas.eps(:,:,1))*1e6,flat(atlas.eps(:,:,2))*1e6,flat(atlas.Gamma_x),flat(atlas.Gamma_y), ...
        flat(atlas.M),1+flat(atlas.Tx),1+flat(atlas.Ty),flat(atlas.R_degree),flat(atlas.R_normalized), ...
        flat(atlas.PRI),flat(atlas.mask99),flat(atlas.pareto), ...
        'VariableNames',{'Phi_deg','Psi_deg','Ex_GPa','Ey_GPa','Gel_GPa','Gx_C_per_m','Gy_C_per_m', ...
        'epsx_400V_microstrain','epsy_400V_microstrain','Gamma_x','Gamma_y','M','A11OC_A11','A22OC_A22', ...
        'R_per_degree','R_normalized','PRI','mask99','pareto'});
    writetable(tableData,fullfile(data,[name '_atlas.csv']));
end
save(fullfile(data,'mode_comparison.mat'),'-struct','modeComparison','-v7');
save(fullfile(data,'mfc_results.mat'),'results','material','summary','extensionSummary', ...
    'benchmarks','losses','convergence','extensionBenchmarks','modeComparison','-v7');
writetable(struct2table(benchmarks),fullfile(data,'benchmark_fingerprints.csv'));
writetable(struct2table(losses),fullfile(data,'finite_family_angle_losses.csv'));
writetable(struct2table(convergence),fullfile(data,'exact_moment_convergence.csv'));
writetable(struct2table(extensionBenchmarks),fullfile(data,'design_extension_benchmarks.csv'));
writetable(angleRows,fullfile(data,'pareto_angle_objectives.csv'));
writetable(frontRows,fullfile(data,'pareto_sum_front.csv'));
mfc_write_json(fullfile(data,'numerical_summary.json'),summary);
mfc_write_json(fullfile(data,'design_extension_summary.json'),extensionSummary);
if makePlots,mfc_plot_all(root);end
record=struct('date',char(datetime('today','Format','yyyy-MM-dd')),'matlab_version',version,'computer',computer, ...
    'elapsed_seconds',toc(started),'entry_point',mfilename,'makePlots',makePlots, ...
    'calculation_reads_python_results',false,'angle_step_deg',.5,'voltage_V',400, ...
    'Vf',.5,'ply_thickness_m',.125e-3,'total_thickness_m',1e-3, ...
    'toolboxes','MATLAB and Optimization Toolbox (fmincon)', ...
    'optimiser','Bounded SQP from six starts; exact 45/45 centre chosen only if tied within 1e-12');
mfc_write_json(fullfile(root,'verification','matlab_run_record.json'),record);
fprintf('Completed MATLAB calculation in %.1f seconds.\n',toc(started));
end
function v=flat(v)
v=reshape(v.',[],1);
end
