function [mask,check] = mfc_pareto(objectives)
%MFC_PARETO Maximise all objectives, preserving equivalent sampled designs.
%   OBJECTIVES has size [ny,nx,nObjectives]. For three objectives the
%   Python tolerance and 12-decimal symmetry equivalence are retained.
shape=size(objectives);flat=reshape(objectives,[],shape(3));
if shape(3)==2
    x=flat(:,1);y=flat(:,2);[~,order]=sortrows([-x,-y],[1 2]);
    keep=false(size(x));ymax=-inf;xtol=1e-10*max(x);ytol=1e-10*max(y);j=1;
    while j<=numel(order)
        k=j+1;
        while k<=numel(order) && abs(x(order(k))-x(order(j)))<=xtol,k=k+1;end
        group=order(j:k-1);mx=max(y(group));
        if mx>ymax+ytol,keep(group(abs(y(group)-mx)<=ytol))=true;end
        ymax=max(ymax,mx);j=k;
    end
    check=struct('nondominated_angle_pair_count',nnz(keep));
else
    tolerance=1e-10;scaled=flat./max(flat,[],1);
    [~,indices,inverse]=unique(round(scaled,12),'rows','first');pool=scaled(indices,:);
    keep=true(size(pool,1),1);
    for start=1:64:size(pool,1)
        ix=start:min(start+63,size(pool,1));
        weak=true(numel(ix),size(pool,1));strict=false(size(weak));
        for k=1:shape(3)
            weak=weak & pool(:,k).'>=pool(ix,k)-tolerance;
            strict=strict | pool(:,k).'>pool(ix,k)+tolerance;
        end
        keep(ix)=~any(weak & strict,2);
    end
    check=struct('normalised_dominance_tolerance',tolerance, ...
        'objective_equivalence_rounding_decimals',12,'unique_objective_count',size(pool,1), ...
        'nondominated_unique_objective_count',nnz(keep), ...
        'nondominated_angle_pair_count',nnz(keep(inverse)));
    % Verify every retained candidate directly against the unrounded grid.
    for candidate=find(keep).'
        dominates=all(scaled>=pool(candidate,:)-tolerance,2) & any(scaled>pool(candidate,:)+tolerance,2);
        assert(~any(dominates),'MFC:Pareto','Retained candidate is dominated.');
    end
    keep=keep(inverse);
end
mask=reshape(keep,shape(1:2));
end
