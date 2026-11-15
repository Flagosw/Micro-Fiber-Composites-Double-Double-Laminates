function r = mfc_score(r,reference)
%MFC_SCORE Apply the unchanged same-mode QI normalisations.
r.Gamma_x=abs(r.G(:,:,1))/abs(reference.G(1));
r.Gamma_y=abs(r.G(:,:,2))/abs(reference.G(2));
r.Sx=r.Ex/reference.Ex;r.Sy=r.Ey/reference.Ey;r.Sgel=r.Gel/reference.Gel;
r.M=hypot(r.Gamma_x./r.Sx,r.Gamma_y./r.Sy);
end
