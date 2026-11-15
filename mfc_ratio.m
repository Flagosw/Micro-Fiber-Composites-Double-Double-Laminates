function [value,guard] = mfc_ratio(numerator,denominator)
%MFC_RATIO Magnitude ratio with the original roundoff-scaled zero guard.
threshold=1e-12*max(max(abs(denominator),[],'all'),1e-300);
valid=abs(denominator)>threshold;value=nan(size(numerator));
value(valid)=abs(numerator(valid))./abs(denominator(valid));
guard=struct('denominator_zero_guard',threshold,'undefined_grid_nodes',nnz(~valid), ...
    'minimum_sampled_abs_denominator',min(abs(denominator),[],'all'), ...
    'signed_denominator_crosses_zero',min(denominator,[],'all')<0 && max(denominator,[],'all')>0, ...
    'signed_numerator_crosses_zero',min(numerator,[],'all')<0 && max(numerator,[],'all')>0);
end
