function enlarge_fonts()
root=fileparts(mfilename('fullpath'));
names={'convergence_mfc','coupling_mfc','free_strain_mfc',...
    'anisotropy_mfc','multifunctionality_mfc','pareto_mfc','pareto_sums_mfc',...
    'mode_comparison_mfc','curvature_mfc','elastic_moduli_mfc',...
    'oc_stiffness_mfc','sensitivity_mfc','pri_mfc'};
out=fullfile(root,'figures'); if ~isfolder(out),mkdir(out);end
report=struct('matlab',version,'font_scale',1.4,'figures',[]);
for k=1:numel(names)
    name=names{k};ref=[name '.pdf'];
    figsrc=fullfile(root,'figures',[name '.fig']);
    f=openfig(figsrc,'invisible'); cleanup=onCleanup(@()close(f));
    ax=findall(f,'Type','axes');
    for a=1:numel(ax),ax(a).XLimMode='manual';ax(a).YLimMode='manual';ax(a).CLimMode='manual';end
    h=findall(f); data={}; fonts={};
    for i=1:numel(h)
        for p={'XData','YData','ZData','CData','CLim','XLim','YLim'}
            if isprop(h(i),p{1}),data(end+1,:)={h(i),p{1},get(h(i),p{1})};end
        end
        if isprop(h(i),'FontSize'),fonts(end+1,:)={h(i),h(i).FontSize};end
    end
    for i=1:size(fonts,1),fonts{i,1}.FontSize=fonts{i,2}*1.4;end
    drawnow;
    for i=1:size(data,1),assert(isequaln(get(data{i,1},data{i,2}),data{i,3}),'Data/limits changed');end
    savefig(f,fullfile(out,[name '.fig']));
    exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white','Padding',8);
    exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180,'BackgroundColor','white','Padding',8);
    r=struct('reference',ref,'source',figsrc,'data_properties_unchanged',size(data,1),'font_objects',size(fonts,1));
    report.figures=[report.figures r];
    fprintf('%s: %d data/limit checks passed\n',name,size(data,1));
    clear cleanup
end
fid=fopen(fullfile(root,'verification','regenerated_font_checks.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
end
