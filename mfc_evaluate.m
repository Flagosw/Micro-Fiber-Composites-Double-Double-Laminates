function r = mfc_evaluate(model,phi,psi,stack,qmean,gmean)
%MFC_EVALUATE Native MATLAB membrane response of the current MFC model.
%   R=MFC_EVALUATE(M,PHI,PSI) uses balanced symmetric DD families.
%   R=MFC_EVALUATE(M,[],[],STACK) uses equal-thickness stack angles.
%   R=MFC_EVALUATE(M,[],[],[],QMEAN,GMEAN) uses prescribed mean tensors.
if nargin>=4 && ~isempty(stack)
    [qb,gb]=mfc_rotate(model.Q,model.g,stack(:).');
    qmean=mean(qb,2);gmean=mean(gb,2);
elseif nargin<5
    [qp,gp]=mfc_rotate(model.Q,model.g,phi);[qm,gm]=mfc_rotate(model.Q,model.g,-phi);
    [rp,hp]=mfc_rotate(model.Q,model.g,psi);[rm,hm]=mfc_rotate(model.Q,model.g,-psi);
    qmean=(qp+qm+rp+rm)/4;gmean=(gp+gm+hp+hm)/4;
else
    qmean=reshape(qmean,1,1,3,3);gmean=reshape(gmean,1,1,3);
end
h=model.h;A=h*qmean;G=h*gmean;K=h*model.kappa;
S=permute(pageinv(permute(A,[3 4 1 2])),[3 4 1 2]);
dA=zeros(size(A));
for i=1:3
    for j=1:3, dA(:,:,i,j)=G(:,:,i).*G(:,:,j)/K; end
end
Aoc=A+dA;Soc=permute(pageinv(permute(Aoc,[3 4 1 2])),[3 4 1 2]);
epsFree=zeros(size(G));
for i=1:3
    for j=1:3,epsFree(:,:,i)=epsFree(:,:,i)+S(:,:,i,j).*G(:,:,j);end
end
efield=model.voltage/model.ell;epsFree=epsFree*efield;
eta=cat(3,G(:,:,1).^2./A(:,:,1,1),G(:,:,2).^2./A(:,:,2,2));
r=struct('A',A,'G',G,'K',K,'Aoc',Aoc,'dA',dA,'S',S,'eps',epsFree, ...
    'Ex',1./(h*S(:,:,1,1)),'Ey',1./(h*S(:,:,2,2)),'Gel',1./(h*S(:,:,3,3)), ...
    'ExOC',1./(h*Soc(:,:,1,1)),'EyOC',1./(h*Soc(:,:,2,2)), ...
    'eta',eta,'Tx',eta(:,:,1)/K,'Ty',eta(:,:,2)/K, ...
    'Ux',.5*eta(:,:,1)*efield^2,'Uy',.5*eta(:,:,2)*efield^2, ...
    'Ufree',.5*efield*sum(G.*epsFree,3));
end
