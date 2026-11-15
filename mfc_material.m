function [models, source] = mfc_material(vf)
%MFC_MATERIAL Calculate active-layer and reduced ply inputs in SI units.
%   [MODELS,SOURCE] = MFC_MATERIAL(VF) implements the mixed-field active
%   layer law from the preserved Python source. Default VF is 0.5.
%   Source: Prasath and Arockiarajan (2014), DOI 10.1016/j.sna.2014.04.008,
%   Table 6, E3 branch. Engineering shear and printed constants are retained.
if nargin==0, vf=0.5; end
sf=compliance([54.05 54.05 48.30]*1e9,[19.48 19.48 19.14]*1e9,[.41 .44 .44]);
sm=compliance([2.9 2.9 2.9]*1e9,[1.11 1.11 1.11]*1e9,[.3 .3 .3]);
df=[-185;-185;440;0;0;0]*1e-12; dm=zeros(6,1);
hf=[sf df;df.' 16.3802e-9]; hm=[sm dm;dm.' .0376e-9];
a=[1 3 5];
k=vf*mixed(hf,a)+(1-vf)*mixed(hm,a);
h=mixed(k,a);
s=h(1:6,1:6); d=h(1:6,7); kt=h(7,7); c=inv(s);
source=struct('Vf_active_layer',vf,'sE_paper_per_Pa',s,'d_E3_paper_m_per_V',d, ...
    'kappaT33_F_per_m',kt,'cE_paper_Pa',c,'e_E3_paper_C_per_m2',c*d, ...
    'kappaS33_F_per_m',kt-d.'*c*d,'full_mixed_law',h, ...
    'constituent_fibre_law',hf,'constituent_matrix_law',hm);
names={'d33','d31'}; indices={[3 2 4],[1 2 6]}; lengths=[380e-6/(1-.2649),125e-6];
for j=1:2
    ix=indices{j}; sr=s(ix,ix); dr=d(ix); q=inv(sr); g=q*dr;
    kp=kt-dr.'*q*dr; qd=q+g*g.'/kp; sd=inv(qd);
    p=struct('paper_voigt_indices_1based',ix,'sE_per_Pa',sr,'QE_Pa',q, ...
        'd_m_per_V',dr,'e_C_per_m2',g,'kappa_plane_strain_F_per_m',kp, ...
        'kappa_T_F_per_m',kt,'QD_Pa',qd,'sD_per_Pa',sd, ...
        'E1E_Pa',1/sr(1,1),'E2E_Pa',1/sr(2,2),'G12E_Pa',1/sr(3,3), ...
        'nu12E',-sr(1,2)/sr(1,1),'nu21E',-sr(1,2)/sr(2,2), ...
        'E1D_Pa',1/sd(1,1),'E2D_Pa',1/sd(2,2),'G12D_Pa',1/sd(3,3), ...
        'nu12D',-sd(1,2)/sd(1,1),'nu21D',-sd(1,2)/sd(2,2));
    source.(names{j})=p;
    models.(names{j})=struct('Q',q,'d',dr,'g',g,'kappa',kp, ...
        'ell',lengths(j),'h',1e-3,'voltage',400,'name',names{j});
end
end

function s=compliance(e,g,nu)
s=diag([1./e,1./g]);
s(1,2)=-nu(1)/e(1);s(2,1)=s(1,2);
s(1,3)=-nu(2)/e(1);s(3,1)=s(1,3);
s(2,3)=-nu(3)/e(2);s(3,2)=s(2,3);
end

function k=mixed(h,a)
r=eye(7);r(a,:)=h(a,:); ri=inv(r);
k=h*ri;k(a,:)=ri(a,:);
end
