function mfc_plot_all(root)
%MFC_PLOT_ALL Render publication figures from native MATLAB calculations.
%   Every plot is exported as editable FIG, PDF and PNG. Python numerical
%   results are never read. Colormaps preserve the original plotting style.
if nargin==0,root=fileparts(mfilename('fullpath'));end
loaded=load(fullfile(root,'data','mfc_results.mat'));
results=loaded.results;summary=loaded.summary;extensionSummary=loaded.extensionSummary;
palette=load(fullfile(root,'colormaps.mat'));
out=fullfile(root,'figures');if ~isfolder(out),mkdir(out);end
names={'d33','d31'};modeLabels={'$d_{33}$','$d_{31}$'};
colors=[0 .447 .698;.835 .369 0;0 .620 .451;.800 .475 .655];

for output={'coupling','oc_stiffness','free_strain'}
    name=output{1};[f,layout]=newFigure(2,2,6.5,5.1);
    for i=1:2
        z=results.(names{i});r=z.atlas;
        for j=1:2
            ax=nexttile(layout);direction=char('x'+j-1);map=palette.viridis;limits=[];
            switch name
                case 'coupling'
                    keys={'Gamma_x','Gamma_y'};v=r.(keys{j});label=['$\Gamma_' direction '$'];
                case 'oc_stiffness'
                    keys={'Tx','Ty'};v=1+r.(keys{j});label=sprintf('$A_{%d%d}^{OC}/A_{%d%d}$',j,j,j,j);
                case 'free_strain'
                    v=r.eps(:,:,j)*1e6;label=['$\varepsilon_' direction '^{free}$ ($\mu\varepsilon$)'];
                    map=palette.Blues_r;
                    if i==1,limits=[-1 1]*max(abs(v),[],'all');map=palette.RdBu_r;end
            end
            mapPanel(ax,r.phi,r.psi,v,sprintf('(%c) %s, %s direction',char(96+2*(i-1)+j),modeLabels{i},direction),label,map,limits,true);
        end
    end
    finish(f,name,out);
end

spec={'multifunctionality','M','$M$';'sensitivity','R_normalized','$\widehat R_M$';'pri','PRI','$PRI$'};
for k=1:size(spec,1)
    [f,layout]=newFigure(1,2,6.5,2.7);
    for i=1:2
        r=results.(names{i}).atlas;
        mapPanel(nexttile(layout),r.phi,r.psi,r.(spec{k,2}),sprintf('(%c) %s',char(96+i),modeLabels{i}),spec{k,3},palette.viridis,[],true);
    end
    finish(f,spec{k,1},out);
end

[f,layout]=newFigure(3,2,6.5,7.5);keys={'Sx','Sy','Sgel'};
labels={'$E_x/E_x^{QI}$','$E_y/E_y^{QI}$','$G_{xy}^{el}/G_{xy}^{el,QI}$'};
for i=1:2
    r=results.(names{i}).atlas;
    for j=1:3,mapPanel(nexttile(layout,2*(j-1)+i),r.phi,r.psi,r.(keys{j}),sprintf('(%c) %s',char(96+2*(j-1)+i),modeLabels{i}),labels{j},palette.viridis,[],true);end
end
finish(f,'elastic_moduli',out);

[f,layout]=newFigure(3,2,6.5,7.5);keys={'A_E','A_G','A_epsilon'};
labels={'$\log_{10}A_E$','$\log_{10}A_G$','$\log_{10}A_\varepsilon$'};
for i=1:2
    r=results.(names{i}).atlas;e=results.(names{i}).extensions;
    for j=1:3
        ax=nexttile(layout,2*(j-1)+i);value=e.(keys{j});value(value<=0)=NaN;value=log10(value);
        peak=max(abs(value),[],'all','omitnan');bound=niceCeiling(peak);
        mapPanel(ax,r.phi,r.psi,value,sprintf('(%c) %s',char(96+2*(j-1)+i),modeLabels{i}),labels{j},palette.RdBu_r,[-bound bound],true);
        if j>1
            if j==2,raw=r.G;else,raw=r.eps;end
            for direction=1:2
                v=raw(:,:,direction);
                if min(v,[],'all')<0 && max(v,[],'all')>0,contour(ax,r.phi,r.psi,v,[0 0],'w--','LineWidth',.8);end
            end
        end
    end
end
finish(f,'anisotropy',out);

[f,layout]=newFigure(1,2,6.5,2.7);mc=loaded.modeComparison;
mapPanel(nexttile(layout),mc.phi,mc.psi,mc.R_G,'(a) x-direction coupling','$R_G$',palette.viridis,[],true);
mapPanel(nexttile(layout),mc.phi,mc.psi,mc.R_epsilon,'(b) x-direction free strain, 400 V','$R_\varepsilon$',palette.viridis,[],true);
finish(f,'mode_comparison',out);
[f,layout]=newFigure(1,2,6.5,2.7);
for i=1:2
    r=results.(names{i}).atlas;e=results.(names{i}).extensions;
    mapPanel(nexttile(layout),r.phi,r.psi,1e3*e.H_M_spectral_radius,sprintf('(%c) %s',char(96+i),modeLabels{i}), ...
        '$10^3\|\mathbf H_M\|_2$ (deg$^{-2}$)',palette.viridis,[],true);
end
finish(f,'curvature',out);

[f,layout]=newFigure(1,2,6.5,3.1);
for i=1:2
    ax=nexttile(layout);z=results.(names{i});r=z.atlas;p=r.pareto;
    gx=abs(r.G(:,:,1))*1e3;
    scatter(ax,r.Ex(:)/1e9,gx(:),3,r.M(:),'filled','MarkerFaceAlpha',.55);hold(ax,'on');colormap(ax,palette.viridis);
    h1=scatter(ax,r.Ex(p)/1e9,gx(p),20,colors(2,:),'filled');
    h2=plot(ax,z.optimum.Ex/1e9,abs(z.optimum.G(1))*1e3,'p','MarkerFaceColor',[.94 .89 .26],'MarkerEdgeColor','k','MarkerSize',9);
    h3=plot(ax,z.reference.Ex/1e9,abs(z.reference.G(1))*1e3,'d','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',5);
    xlabel(ax,'$E_x$ (GPa)','Interpreter','latex');ylabel(ax,'$|G_x|$ ($10^{-3}$ C/m)','Interpreter','latex');
    title(ax,sprintf('(%c) %s',char(96+i),modeLabels{i}),'Interpreter','latex');
    cb=colorbar(ax);cb.Label.String='$M$';cb.Label.Interpreter='latex';cb.TickLabelInterpreter='latex';
    legend(ax,[h1 h2 h3],{'Nondominated','$M$ maximum','QI = QUAD'},'Interpreter','latex','FontSize',6.5,'Location','best','Box','off');
    styleAxes(ax);grid(ax,'off');ax.GridAlpha=.12;
end
finish(f,'pareto',out);

[f,layout]=newFigure(2,2,7.2,5.8);
for i=1:2
    z=results.(names{i});r=z.atlas;e=z.extensions;p=e.pareto3;ax=nexttile(layout,i);
    xx=e.E_sum/1e9;yy=e.strain_abs_sum*1e6;
    if extensionSummary.(names{i}).G_abs_sum_relative_span<1e-10
        scatter(ax,xx(:),yy(:),3,[.50 .60 .67],'filled','MarkerFaceAlpha',.45);
    else
        gref=sum(abs(z.reference.G(1:2)));cc=e.G_abs_sum/gref;
        scatter(ax,xx(:),yy(:),3,cc(:),'filled','MarkerFaceAlpha',.55);colormap(ax,palette.viridis);
        cb=colorbar(ax);cb.Label.String='$(|G_x|+|G_y|)/G_\Sigma^{QI}$';cb.Label.Interpreter='latex';cb.TickLabelInterpreter='latex';
    end
    hold(ax,'on');h=scatter(ax,xx(p),yy(p),10,colors(2,:),'filled');
    points={z.reference,z.dd45,z.ud};markers={'d','p','o'};pointColors={'w',[.94 .89 .26],colors(3,:)};pointNames={'QI','DD45','UDx'};
    for k=1:3
        b=points{k};px=(b.Ex+b.Ey)/1e9;py=sum(abs(b.eps(1:2)))*1e6;
        plot(ax,px,py,markers{k},'MarkerSize',6,'MarkerFaceColor',pointColors{k},'MarkerEdgeColor','k');
        text(ax,px,py,[' ' pointNames{k}],'FontSize',6.5,'VerticalAlignment','bottom');
    end
    xlabel(ax,'$E_x+E_y$ (GPa)','Interpreter','latex');ylabel(ax,'$|\varepsilon_x^{free}|+|\varepsilon_y^{free}|$ ($\mu\varepsilon$)','Interpreter','latex');
    title(ax,sprintf('(%c) %s: objective projection',char(96+i),modeLabels{i}),'Interpreter','latex');
    legend(ax,h,'Nondominated (3 objectives)','FontSize',6.5,'Location','southwest','Box','off');styleAxes(ax);grid(ax,'off');ax.GridAlpha=.12;
    ax=nexttile(layout,i+2);scatter(ax,r.phi(:),r.psi(:),1,[.88 .88 .88],'filled');hold(ax,'on');
    scatter(ax,r.phi(p),r.psi(p),18,colors(2,:),'filled','Clipping','off');marks(ax,true);axis(ax,'square');
    title(ax,sprintf('(%c) %s: nondominated angle pairs',char(98+i),modeLabels{i}),'Interpreter','latex');
    text(ax,3,86,'$(0^\circ,90^\circ)$','Interpreter','latex','FontSize',7);
    text(ax,87,4,'$(90^\circ,0^\circ)$','Interpreter','latex','FontSize',7,'HorizontalAlignment','right');
end
finish(f,'pareto_sums',out);

[f,layout]=newFigure(2,2,6.5,5.2);
for i=1:2
    z=results.(names{i});[xp,yp]=meshgrid(linspace(-10,10,201));
    val=mfc_score(mfc_evaluate(z.model,z.opt(1)+xp,z.opt(2)+yp),z.reference);
    loss=100*(1-val.M/summary.(names{i}).Mmax);ax=nexttile(layout);
    mapPanel(ax,xp,yp,loss,sprintf('(%c) %s: local loss',char(95+2*i),modeLabels{i}), ...
        '$100(1-M/M_{max})$ (\%)',palette.magma_r,[],false);
    xlim(ax,[-10 10]);ylim(ax,[-10 10]);xticks(ax,[-10 -5 0 5 10]);yticks(ax,[-10 -5 0 5 10]);
    xlabel(ax,'$\Delta\Phi$ (deg)','Interpreter','latex');ylabel(ax,'$\Delta\Psi$ (deg)','Interpreter','latex');
    plot(ax,0,0,'p','MarkerFaceColor',[.94 .89 .26],'MarkerEdgeColor','k','MarkerSize',7);
    ax=nexttile(layout);rr=summary.(names{i}).finite_family_angle_losses;
    plot(ax,[rr.delta_deg],100*[rr.worst_loss_fraction],'o-','Color',colors(2,:),'LineWidth',1.2);hold(ax,'on');
    plot(ax,[rr.delta_deg],100*[rr.mean_loss_fraction],'s-','Color',colors(1,:),'LineWidth',1.2);
    xlabel(ax,'Whole-family bound $\delta$ (deg)','Interpreter','latex');ylabel(ax,'Score loss (%)');
    title(ax,sprintf('(%c) %s: finite errors',char(96+2*i),modeLabels{i}),'Interpreter','latex');
    xticks(ax,[1 2 5 10]);legend(ax,{'Worst sampled loss','Mean sampled loss'},'Location','northwest','FontSize',7,'Box','off');styleAxes(ax);grid(ax,'off');
end
finish(f,'local_robustness',out);

[f,layout]=newFigure(2,2,6.5,5.1);intervals=[20 70;10 50];
keys={'error_A','error_G_vector','error_H400','error_U_OC'};
labels={'$\eta_A$','$\eta_G$','$\eta_{\mathcal H(400V)}$','$\eta_{U_{OC}}$'};styles={'-','--','-.',':'};
for i=1:2
    for j=1:2
        ax=nexttile(layout);allRows=loaded.convergence;
        cr=allRows(strcmp({allRows.mode},names{i}) & [allRows.low_angle_deg]==intervals(j,1));
        curves=gobjects(4,1);
        for k=1:4,curves(k)=semilogy(ax,[cr.N_positive],max([cr.(keys{k})],1e-16),styles{k},'Color',colors(k,:),'LineWidth',1.2);hold(ax,'on');end
        yline(ax,.01,':','Color',[.5 .5 .5],'HandleVisibility','off');
        if j==1,ylim(ax,[1e-16 1]);else,ylim(ax,[1e-5 1]);end
        xlim(ax,[2 200]);xlabel(ax,'Positive angular states $N$','Interpreter','latex');ylabel(ax,'Relative error');
        title(ax,sprintf('(%c) %s: $%d^\\circ$--$%d^\\circ$',char(96+2*(i-1)+j),modeLabels{i},intervals(j,:)),'Interpreter','latex');
        location='northeast';if j==1,location='east';end
        legend(ax,curves,labels,'Interpreter','latex','FontSize',6.5,'NumColumns',2,'Location',location,'Box','off');styleAxes(ax);grid(ax,'off');ax.GridAlpha=.12;
    end
end
finish(f,'convergence',out);

% Additional figures produced by the Python analysis but not cited in the
% current manuscript are retained locally for full reproduction.
[f,layout]=newFigure(1,2,6.5,2.8);
for i=1:2
    r=results.(names{i}).atlas;rel=r.M/summary.(names{i}).Mmax;rel(~r.mask99)=NaN;
    ax=nexttile(layout);mapPanel(ax,r.phi,r.psi,rel,sprintf('(%c) %s',char(96+i),modeLabels{i}), ...
        '$M/M_{max}\geq0.99$',palette.viridis,[.99 1],true);ax.Color=[.94 .94 .94];
end
finish(f,'near_optimal',out);
[f,layout]=newFigure(1,2,6.5,3.2);
metrics={'Sx','Gamma_x','abs_epsx_ratio','Tx_ratio','Ux_ratio'};
for i=1:2
    ax=nexttile(layout);b=loaded.benchmarks;values=zeros(5,3);arch={'QI','QUAD','DD optimum'};
    for j=1:3
        row=b(strcmp({b.mode},names{i}) & strcmp({b.architecture},arch{j}));
        for k=1:5,values(k,j)=row.(metrics{k});end
    end
    bars=bar(ax,values,'grouped');for j=1:3,bars(j).FaceColor=colors(j,:);end
    ax.XTickLabel={'$S_x$','$\Gamma_x$','$|\varepsilon_x|$ ratio','$\bar T_x$ ratio','$U_e^{(x)}$ ratio'};
    ax.XTickLabelRotation=25;ax.TickLabelInterpreter='latex';ylim(ax,[0 1.55]);
    ylabel(ax,'Ratio to same-mode QI');title(ax,sprintf('(%c) %s',char(96+i),modeLabels{i}),'Interpreter','latex');
    legend(ax,arch,'Location','northwest','FontSize',6.5,'Box','off');styleAxes(ax);yline(ax,1,'k:','HandleVisibility','off');
end
finish(f,'fingerprint',out);
% Figure 1 is separately authored schematic artwork.
fprintf('Exported 16 quantitative figure sets.\n');
end

function [f,layout]=newFigure(rows,cols,width,height)
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 width height]);
layout=tiledlayout(f,rows,cols,'TileSpacing','compact','Padding','compact');
end
function styleAxes(ax)
set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.6,'TickLabelInterpreter','latex','Box','on','Layer','top');
end
function marks(ax,show)
if show
    plot(ax,22.5,67.5,'d','MarkerSize',4,'MarkerFaceColor','w','MarkerEdgeColor','k');
    plot(ax,45,45,'p','MarkerSize',7,'MarkerFaceColor',[.94 .89 .26],'MarkerEdgeColor','k');
end
xlim(ax,[0 90]);ylim(ax,[0 90]);xticks(ax,[0 30 60 90]);yticks(ax,[0 30 60 90]);
xlabel(ax,'$\Phi$ (deg)','Interpreter','latex');ylabel(ax,'$\Psi$ (deg)','Interpreter','latex');styleAxes(ax);
end
function mapPanel(ax,x,y,value,titleText,label,cmap,limits,showMarks)
im=imagesc(ax,x(1,:),y(:,1),value);set(im,'AlphaData',isfinite(value));axis(ax,'xy');hold(ax,'on');
colormap(ax,cmap);if ~isempty(limits),clim(ax,limits);end
finite=value(isfinite(value));
if ~isempty(finite) && max(finite)-min(finite)>1e-12*max([1;abs(finite(:))])
    levels=linspace(min(finite),max(finite),7);levels=round(levels(2:6),4,'significant');levels=unique(levels);
    if isscalar(levels),levels=[levels levels];end
    if max(abs(finite))>2 && contains(label,'log'),levels=[-1 -.5 0 .5 1];end
    [cm,~]=contour(ax,x,y,value,levels,'LineColor',[.17 .17 .17],'LineWidth',.45);
    labelContours(ax,cm,x,y);
end
marks(ax,showMarks);if ~showMarks,xlim(ax,[min(x,[],'all') max(x,[],'all')]);ylim(ax,[min(y,[],'all') max(y,[],'all')]);end
axis(ax,'square');title(ax,titleText,'Interpreter','latex','FontSize',9);
cb=colorbar(ax);cb.Label.String=label;cb.Label.Interpreter='latex';cb.TickLabelInterpreter='latex';cb.FontSize=7;cb.Label.FontSize=8;
end
function labelContours(ax,cm,x,y)
% Put one label per true level at an interior point, away from markers.
xlo=min(x,[],'all');xspan=max(x,[],'all')-xlo;
ylo=min(y,[],'all');yspan=max(y,[],'all')-ylo;
levels=[];segments={};j=1;placed=[];
while j<size(cm,2)
    level=cm(1,j);count=cm(2,j);ix=find(levels==level,1);
    if isempty(ix),levels(end+1)=level;segments{end+1}=[];ix=numel(levels);end %#ok<AGROW>
    segments{ix}=[segments{ix},cm(:,j+1:j+count)];j=j+count+1;
end
for k=1:numel(levels)
    points=segments{k};scaled=(points-[xlo;ylo])./[xspan;yspan];
    valid=all(scaled>.1 & scaled<.9,1) & vecnorm(scaled-[.5;.5])>.12 & vecnorm(scaled-[.25;.75])>.14;
    for previous=1:size(placed,2),valid=valid & vecnorm(scaled-placed(:,previous))>.16;end
    points=points(:,valid);scaled=scaled(:,valid);if isempty(points),continue;end
    target=[.28+.44*mod(k,2);.28+.2*mod(k,3)];[~,which]=min(vecnorm(scaled-target));
    placed(:,end+1)=scaled(:,which); %#ok<AGROW>
    label=sprintf('%.3g',levels(k));if abs(levels(k))<1e-12*max(abs(levels)),label='0';end
    cmap=colormap(ax); crange=clim(ax); ci=max(1,min(size(cmap,1),1+round((levels(k)-crange(1))/diff(crange)*(size(cmap,1)-1))));
    textColor=[.05 .05 .05]; if dot(cmap(ci,:),[.2126 .7152 .0722])<.45,textColor=[1 1 1];end
    text(ax,points(1,which),points(2,which),label, ...
        'FontName','Times New Roman','FontSize',6.5,'FontWeight','bold','Color',textColor, ...
        'BackgroundColor','none','Margin',.4,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on');
end
end
function value=niceCeiling(value)
scale=10^floor(log10(max(value,1e-12)));value=ceil(value/scale*2)*scale/2;
end
function finish(f,name,out)
clean=onCleanup(@()close(f));drawnow;
set(findall(f,'Type','axes'),'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off');
savefig(f,fullfile(out,[name '_mfc.fig']));
exportgraphics(f,fullfile(out,[name '_mfc.pdf']),'ContentType','vector','BackgroundColor','white','Padding',5);
exportgraphics(f,fullfile(out,[name '_mfc.png']),'Resolution',220,'BackgroundColor','white','Padding',5);
end
