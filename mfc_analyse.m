function [results,summaries,benchmarks,losses,convergence] = mfc_analyse(models)
%MFC_ANALYSE Reproduce the complete baseline MFC analyses natively.
angle=0:.5:90;[phi,psi]=meshgrid(angle,angle);
names={'d33','d31'};benchmarks=struct([]);losses=struct([]);convergence=struct([]);
for index=1:2
    name=names{index};model=models.(name);
    reference=mfc_evaluate(model,22.5,67.5);reference=mfc_score(reference,reference);
    r=mfc_score(mfc_evaluate(model,phi,psi),reference);
    [~,imax]=max(r.M,[],'all','linear');
    starts=[phi(imax) psi(imax);45 45;0 0;90 90;0 90;22.5 67.5];
    opts=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
        'OptimalityTolerance',1e-12,'StepTolerance',1e-10, ...
        'FiniteDifferenceType','central','MaxFunctionEvaluations',3000);
    solutions=zeros(6,5);
    for k=1:size(starts,1)
        objective=@(x)scoreObjective(model,reference,x);
        assert(isfinite(objective(starts(k,:))));
        [x,value,flag,output]=fmincon(objective,starts(k,:),[],[],[],[],[0 0],[90 90],[],opts);
        assert(flag>0 && all(x>=0 & x<=90),'MFC:Optimisation','Local search failed.');
        solutions(k,:)=[x,value,flag,output.firstorderopt];
    end
    [best,which]=min(solutions(:,3));opt=solutions(which,1:2);
    r45=mfc_score(mfc_evaluate(model,45,45),reference);
    % Pick the exact symmetry centre only when it ties the best local value.
    % This removes insignificant platform-dependent optimum-coordinate noise.
    if abs(best+r45.M)<1e-12,opt=[45 45];end
    ropt=mfc_score(mfc_evaluate(model,opt(1),opt(2)),reference);maxval=ropt.M;
    assert(maxval>=max(r.M,[],'all')-1e-12,'MFC:Optimisation','Search below grid maximum.');
    quad=mfc_score(mfc_evaluate(model,[],[],[0 45 -45 90 90 -45 45 0]),reference);
    ud=mfc_score(mfc_evaluate(model,0,0),reference);
    [dx,dy]=gradientSecondOrder(r.M,.5);
    r.R_degree=hypot(dx,dy);r.R_normalized=r.R_degree/max(r.R_degree,[],'all');
    r.PRI=r.M./(1+r.R_normalized);r.mask99=r.M>=.99*maxval;
    [r.pareto,~]=mfc_pareto(cat(3,r.Ex,abs(r.G(:,:,1))));
    labels={'QI','QUAD','DD45','DD optimum','UD'};values={reference,quad,r45,ropt,ud};
    for k=1:5
        b=values{k};row=struct('mode',name,'architecture',labels{k}, ...
            'Ex_GPa',b.Ex/1e9,'Ey_GPa',b.Ey/1e9,'Gel_GPa',b.Gel/1e9, ...
            'Gx_C_per_m',b.G(1),'Gy_C_per_m',b.G(2), ...
            'epsx_400V_microstrain',b.eps(1)*1e6,'epsy_400V_microstrain',b.eps(2)*1e6, ...
            'Sx',b.Sx,'Gamma_x',b.Gamma_x,'abs_epsx_ratio',abs(b.eps(1)/reference.eps(1)), ...
            'Tx_ratio',b.Tx/reference.Tx,'Ux_ratio',b.Ux/reference.Ux,'M',b.M, ...
            'Tx',b.Tx,'Ty',b.Ty,'Ux_J_per_m2',b.Ux,'Ufree_J_per_m2',b.Ufree, ...
            'ExOC_GPa',b.ExOC/1e9,'EyOC_GPa',b.EyOC/1e9);
        benchmarks=[benchmarks;row]; %#ok<AGROW>
    end
    modeLosses=struct([]);
    for delta=[1 2 5 10]
        dd=linspace(-delta,delta,201);[dp,ds]=meshgrid(dd,dd);
        value=mfc_score(mfc_evaluate(model,opt(1)+dp,opt(2)+ds),reference);
        val=value.M;
        row=struct('mode',name,'delta_deg',delta,'nominal_phi_deg',opt(1), ...
            'nominal_psi_deg',opt(2),'worst_loss_fraction',1-min(val,[],'all')/maxval, ...
            'mean_loss_fraction',1-mean(val,'all')/maxval, ...
            'population_CV',std(val(:),1)/mean(val,'all'),'sample_count',numel(val));
        modeLosses=[modeLosses;row];losses=[losses;row]; %#ok<AGROW>
    end
    aScale=max(abs(reference.A),[],'all');
    eigenA=matrixEigenvalues(r.A);eigenIncrement=matrixEigenvalues(r.dA);
    validations=struct('QI_QUAD_A_relative_error',max(abs(reference.A-quad.A),[],'all')/aScale, ...
        'QI_QUAD_G_relative_error',max(abs(reference.G-quad.G),[],'all')/max(abs(reference.G),[],'all'), ...
        'QI_inplane_isotropy_relative_error',max(abs(reference.A(1,1,1,1)-reference.A(1,1,2,2)), ...
        abs(reference.A(1,1,3,3)-(reference.A(1,1,1,1)-reference.A(1,1,1,2))/2))/aScale, ...
        'minimum_A_eigenvalue_N_per_m',min(eigenA,[],'all'), ...
        'minimum_OC_increment_eigenvalue_N_per_m',min(eigenIncrement,[],'all'), ...
        'max_balance_residual_N_per_m',max(abs(r.A(:,:,1:2,3)),[],'all'), ...
        'max_Gxy_C_per_m',max(abs(r.G(:,:,3)),[],'all'), ...
        'phi_psi_M_symmetry_error',max(abs(r.M-r.M.'),[],'all'), ...
        'x_y_complement_M_symmetry_error',max(abs(r.M-flip(flip(r.M,1),2)),[],'all'), ...
        'min_constrained_OC_ratio',min([1+r.Tx(:);1+r.Ty(:)]), ...
        'min_Young_modulus_OC_ratio',min([r.ExOC(:)./r.Ex(:);r.EyOC(:)./r.Ey(:)]), ...
        'covariance_positive',model.kappa>0,'field_400V_V_per_m',model.voltage/model.ell);
    assert(validations.minimum_A_eigenvalue_N_per_m>0);
    assert(validations.QI_QUAD_A_relative_error<1e-12);
    assert(validations.minimum_OC_increment_eigenvalue_N_per_m>-1e-8);
    assert(validations.min_Young_modulus_OC_ratio>=1-1e-12);
    xx=r.Ex;yy=abs(r.G(:,:,1));xtol=1e-10*max(xx,[],'all');ytol=1e-10*max(yy,[],'all');
    s=struct('optimum_phi_deg',opt(1),'optimum_psi_deg',opt(2),'Mmax',maxval, ...
        'M_DD45',r45.M,'M_UD',ud.M, ...
        'near99_phi_marginal_deg',[min(phi(r.mask99)) max(phi(r.mask99))], ...
        'near99_psi_marginal_deg',[min(psi(r.mask99)) max(psi(r.mask99))], ...
        'near99_fraction_grid',mean(r.mask99,'all'),'pareto_grid_count',nnz(r.pareto), ...
        'pareto_phi_psi_deg',[phi(r.pareto) psi(r.pareto)], ...
        'optimum_dominated_by_grid',any(xx>=ropt.Ex-xtol & yy>=abs(ropt.G(1))-ytol & ...
        (xx>ropt.Ex+xtol | yy>abs(ropt.G(1))+ytol),'all'), ...
        'free_strain_x_range_microstrain',[min(r.eps(:,:,1),[],'all') max(r.eps(:,:,1),[],'all')]*1e6, ...
        'free_strain_y_range_microstrain',[min(r.eps(:,:,2),[],'all') max(r.eps(:,:,2),[],'all')]*1e6, ...
        'validations',validations,'finite_family_angle_losses',modeLosses);
    intervals=[20 70;10 50];convSummary=struct([]);
    for j=1:2
        lo=intervals(j,1);hi=intervals(j,2);
        [qe,ge,mom]=mfc_exact_spectrum(model,lo,hi);
        [he,hsc,uoce,scale]=mfc_energy(qe,ge,model.kappa,model.voltage/model.ell);
        intervalRows=struct([]);
        for n=2:2:200
            tt=linspace(lo,hi,n);[qb,gb]=mfc_rotate(model.Q,model.g,[tt -tt]);
            qn=squeeze(mean(qb,2));gn=squeeze(mean(gb,2));
            [hn,hscn,uocn]=mfc_energy(qn,gn,model.kappa,model.voltage/model.ell);
            hError=NaN;if abs(he)>1e-10*scale,hError=abs(hn-he)/abs(he);end
            row=struct('mode',name,'low_angle_deg',lo,'high_angle_deg',hi,'N_positive',n, ...
                'error_A',norm(qn-qe,'fro')/norm(qe,'fro'), ...
                'error_G_vector',norm(gn-ge)/norm(ge),'error_H400',hError, ...
                'error_H0',abs(hscn-hsc)/abs(hsc),'error_U_OC',abs(uocn-uoce)/abs(uoce), ...
                'error_H_positive_scale',abs(hn-he)/scale,'H400_J_per_m3',hn, ...
                'H400_exact_J_per_m3',he,'population_variance_gx',var(gb(:,:,1),1,2));
            intervalRows=[intervalRows;row];convergence=[convergence;row]; %#ok<AGROW>
        end
        keys={'error_A','error_G_vector','error_H400','error_H0','error_U_OC'};thresholds=struct;
        for k=1:numel(keys)
            at=find([intervalRows.(keys{k})]<.01,1);
            if isempty(at),thresholds.(keys{k})=NaN;else,thresholds.(keys{k})=intervalRows(at).N_positive;end
        end
        item=struct('interval_deg',[lo hi],'exact_Q_Pa',qe,'exact_g_C_per_m2',ge, ...
            'exact_moments',mom,'H400_exact_J_per_m3',he,'threshold_1percent',thresholds);
        convSummary=[convSummary;item]; %#ok<AGROW>
    end
    s.convergence=convSummary;summaries.(name)=s;
    r.phi=phi;r.psi=psi;
    results.(name)=struct('model',model,'atlas',r,'reference',reference,'opt',opt, ...
        'optimum',ropt,'quad',quad,'dd45',r45,'ud',ud,'optimisation_runs',solutions);
    fprintf('%s: Mmax %.12g at %.9f/%.9f degrees; Pareto %d\n',name,maxval,opt,nnz(r.pareto));
end
end

function val=scoreObjective(model,reference,x)
r=mfc_score(mfc_evaluate(model,x(1),x(2)),reference);val=-r.M;
end

function [dx,dy]=gradientSecondOrder(z,h)
[dx,dy]=gradient(z,h);
dx(:,1)=(-3*z(:,1)+4*z(:,2)-z(:,3))/(2*h);
dx(:,end)=(3*z(:,end)-4*z(:,end-1)+z(:,end-2))/(2*h);
dy(1,:)=(-3*z(1,:)+4*z(2,:)-z(3,:))/(2*h);
dy(end,:)=(3*z(end,:)-4*z(end-1,:)+z(end-2,:))/(2*h);
end

function values=matrixEigenvalues(a)
pages=reshape(permute(a,[3 4 1 2]),3,3,[]);values=zeros(size(pages,3),3);
for k=1:size(pages,3),values(k,:)=sort(eig((pages(:,:,k)+pages(:,:,k).')/2)).';end
end
