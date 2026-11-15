function h = mfc_hessian(model,reference,phi,psi,step)
%MFC_HESSIAN Central angular differences with derivatives per degree^2.
score=@(p,s)localScore(model,reference,p,s);
center=score(phi,psi);
h=zeros([size(center),2,2]);
h(:,:,1,1)=(score(phi+step,psi)-2*center+score(phi-step,psi))/step^2;
h(:,:,2,2)=(score(phi,psi+step)-2*center+score(phi,psi-step))/step^2;
h(:,:,1,2)=(score(phi+step,psi+step)-score(phi+step,psi-step) ...
    -score(phi-step,psi+step)+score(phi-step,psi-step))/(4*step^2);
h(:,:,2,1)=h(:,:,1,2);
end
function value=localScore(model,reference,phi,psi)
r=mfc_score(mfc_evaluate(model,phi,psi),reference);value=r.M;
end
