function revise_review_figures()
% Revise only figure presentation requested by the 27 October comments.
root=fileparts(mfilename('fullpath'));
base=root;
out=fullfile(root,'figures');
if ~isfolder(out),mkdir(out);end
saved=load(fullfile(base,'data','mfc_results.mat'));
pal=load(fullfile(base,'colormaps.mat'));
names={'d33','d31'};

% Retain the plotted data and axes; replace the classification legend.
f=openfig(fullfile(base,'figures','pareto_mfc.fig'),'invisible');
allLegends=findall(f,'Type','legend');
for k=1:numel(allLegends)
    allLegends(k).String={'$(0^\circ,0^\circ)$','$(45^\circ,45^\circ)$','$\mathrm{QI}\;(22.5^\circ,67.5^\circ)$'};
    allLegends(k).Interpreter='latex';
    allLegends(k).FontSize=6.5;
end
finish(f,'pareto_mfc',out);

% Use the same complete sample arrays, without Pareto overlays or angle maps.
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.2 3.3]);
t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for i=1:2
    z=saved.results.(names{i}); e=z.extensions; ax=nexttile(t);
    xx=e.E_sum/1e9; yy=e.strain_abs_sum*1e6;
    if i==1
        cc=e.G_abs_sum/sum(abs(z.reference.G(1:2)));
        scatter(ax,xx(:),yy(:),3,cc(:),'filled','MarkerFaceAlpha',.55);
        colormap(ax,pal.viridis);
        cb=colorbar(ax);cb.Label.String='$G_\Sigma/G_\Sigma^{QI}$';
        cb.Label.Interpreter='latex';cb.TickLabelInterpreter='latex';cb.FontSize=7;
    else
        scatter(ax,xx(:),yy(:),3,[.50 .60 .67],'filled','MarkerFaceAlpha',.45);
    end
    hold(ax,'on');
    pts={z.reference,z.dd45,z.ud}; markers={'d','p','o'};
    colors={'w',[.94 .89 .26],[0 .62 .451]};
    labels={'$\mathrm{QI}$','$(45^\circ,45^\circ)$','$\mathrm{UD}_x$'};
    xrange=max(xx,[],'all')-min(xx,[],'all');
    yrange=max(yy,[],'all')-min(yy,[],'all');
    for k=1:3
        b=pts{k};px=(b.Ex+b.Ey)/1e9;py=sum(abs(b.eps(1:2)))*1e6;
        plot(ax,px,py,markers{k},'MarkerSize',7,'MarkerFaceColor',colors{k},'MarkerEdgeColor','k');
        tx=px+.025*xrange;ty=py+.025*yrange;ha='left';va='bottom';
        if k==3,tx=px-.025*xrange;ha='right';end
        if i==2 && k<=2,ty=py-.045*yrange;va='top';end
        text(ax,tx,ty,labels{k},'Interpreter','latex','FontSize',7,'HorizontalAlignment',ha,'VerticalAlignment',va);
    end
    xlabel(ax,'$E_\Sigma$ (GPa)','Interpreter','latex');
    ylabel(ax,'$\varepsilon_\Sigma$ ($\mu\varepsilon$)','Interpreter','latex');
    title(ax,sprintf('(%c)',96+i),'Interpreter','latex');
    ylim(ax,[min(yy,[],'all')-.04*yrange,max(yy,[],'all')+.08*yrange]);
    set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.6,...
        'TickLabelInterpreter','latex','Box','on','Layer','top');
    grid(ax,'off');
end
finish(f,'pareto_sums_mfc',out);

% Titles are presentation only; all map values and contour labels are kept.
f=openfig(fullfile(base,'figures','mode_comparison_mfc.fig'),'invisible');
aa=findall(f,'Type','axes');
for k=1:numel(aa)
    s=string(aa(k).Title.String);
    if startsWith(s,'(a)'),aa(k).Title.String='(a)';end
    if startsWith(s,'(b)'),aa(k).Title.String='(b)';end
end
finish(f,'mode_comparison_mfc',out);
fprintf('Exported three revised figures; numerical arrays are unchanged.\n');
end

function finish(f,name,out)
clean=onCleanup(@()close(f));drawnow;
set(findall(f,'Type','axes'),'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off');
savefig(f,fullfile(out,[name '.fig']));
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white','Padding',5);
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',220,'BackgroundColor','white','Padding',5);
end
