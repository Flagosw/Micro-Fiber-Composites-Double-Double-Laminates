function [enthalpy,mechanical,uoc,scale] = mfc_energy(q,g,kappa,efield)
%MFC_ENERGY Electric enthalpy and SC/OC energies, all in J/m^3.
strain=[.001;0;0];mechanical=.5*strain.'*q*strain;
enthalpy=mechanical-strain.'*g*efield-.5*kappa*efield^2;
uoc=.5*strain.'*(q+g*g.'/kappa)*strain;
scale=mechanical+abs(strain.'*g*efield)+.5*kappa*efield^2;
end
