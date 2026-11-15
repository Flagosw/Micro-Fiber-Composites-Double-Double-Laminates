function revise_legends()
% Add explicit angle legends to the current, preserved MFC figures.
root=fileparts(mfilename('fullpath'));
base=root;
recent=fullfile(root,'figures');
out=fullfile(root,'figures');
if ~isfolder(out),mkdir(out);end
oldPath=path;cleanup=onCleanup(@()path(oldPath));
code=root;
addpath(code);
assert(strcmp(which('mfc_evaluate'),fullfile(code,'mfc_evaluate.m')));
report=struct('matlab',version,'source_model',code,'figures',struct());
names={'coupling','free_strain','anisotropy','multifunctionality',...
    'mode_comparison','curvature','elastic_moduli','oc_stiffness','sensitivity','pri'};
for j=1:numel(names)
    name=[names{j} '_mfc'];
    src=fullfile(base,'figures',[name '.fig']);
    if strcmp(names{j},'mode_comparison'),src=fullfile(recent,[name '.fig']);end
    f=openfig(src,'invisible');
    old=recordData(f);
    layout=findall(f,'Type','tiledlayout');
    ax=findall(f,'Type','axes');
    assert(isscalar(layout));
    diamond=findobj(ax(1),'Type','line','Marker','d');
    star=findobj(ax(1),'Type','line','Marker','p');
    assert(isscalar(diamond)&&isscalar(star));
    assert(diamond.XData==22.5 && diamond.YData==67.5);
    assert(star.XData==45 && star.YData==45);
    lg=legend(ax(1),[diamond star],...
        {'$\mathrm{QI}:\ (\Phi,\Psi)=(22.5^\circ,67.5^\circ)$',...
        '$(\Phi,\Psi)=(45^\circ,45^\circ)$'},...
        'Interpreter','latex','FontSize',8,'NumColumns',2,'Box','off');
    lg.Layout.Tile='south';lg.AutoUpdate='off';
    f.Units='inches';p=f.Position;p(4)=p(4)+.32;f.Position=p;
    assertSame(old);
    finish(f,name,out);
    report.figures.(names{j})=struct('source',src,'data_properties_unchanged',size(old,1),...
        'legend_angles_deg',[22.5 67.5;45 45]);
end

% The summed-objective plot has three architectures, including UDx.
src=fullfile(recent,'pareto_sums_mfc.fig');f=openfig(src,'invisible');
old=recordData(f);aa=findall(f,'Type','axes');
for k=1:numel(aa)
    labels=findall(aa(k),'Type','text');
    for n=1:numel(labels)
        if any(strcmp(labels(n).String,{'$\mathrm{QI}$','$(45^\circ,45^\circ)$','$\mathrm{UD}_x$'}))
            labels(n).Visible='off';
        end
    end
end
ax=aa(1);d=findobj(ax,'Type','line','Marker','d');
p=findobj(ax,'Type','line','Marker','p');o=findobj(ax,'Type','line','Marker','o');
lg=legend(ax,[d p o],{'$\mathrm{QI}:\ (22.5^\circ,67.5^\circ)$',...
    '$(45^\circ,45^\circ)$','$\mathrm{UD}_x:\ (0^\circ,0^\circ)$'},...
    'Interpreter','latex','FontSize',8,'NumColumns',3,'Box','off');
lg.Layout.Tile='south';lg.AutoUpdate='off';
f.Units='inches';pos=f.Position;pos(4)=pos(4)+.32;f.Position=pos;
assertSame(old);finish(f,'pareto_sums_mfc',out);
report.figures.pareto_sums=struct('source',src,'data_properties_unchanged',size(old,1));

% Independent closed form and shared-model check of the complementary line.
models=mfc_material(.5);phi=0:.125:90;
for mode={'d33','d31'}
    m=models.(mode{1});r=mfc_evaluate(m,phi,90-phi);
    expected=sum(m.g(1:2))/(m.Q(1,1)+m.Q(2,2)+2*m.Q(1,2))*m.voltage/m.ell;
    residual=max(abs(r.eps(:,:,1:2)-expected),[],'all');
    assert(residual<1e-15);
    report.complementary_line.(mode{1})=struct('samples',numel(phi),...
        'voltage_V',m.voltage,'strain_microstrain',expected*1e6,...
        'closed_form_max_residual',residual,'shear_max',max(abs(r.eps(:,:,3)),[],'all'));
end
mfc_write_json(fullfile(root,'verification','regenerated_legend_checks.json'),report);
disp(report.complementary_line);
fprintf('Eleven figures exported; every original X/Y/Z/C data property unchanged.\n');
end

function record=recordData(f)
record={};all=findall(f);props={'XData','YData','ZData','CData'};
for i=1:numel(all)
    for j=1:numel(props)
        if isprop(all(i),props{j})
            record(end+1,:)={all(i),props{j},get(all(i),props{j})}; %#ok<AGROW>
        end
    end
end
end
function assertSame(record)
for i=1:size(record,1)
    assert(isequaln(get(record{i,1},record{i,2}),record{i,3}),'Plot data changed.');
end
end
function finish(f,name,out)
cleanup=onCleanup(@()close(f));drawnow;
savefig(f,fullfile(out,[name '.fig']));
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white','Padding',5);
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',220,'BackgroundColor','white','Padding',5);
end
