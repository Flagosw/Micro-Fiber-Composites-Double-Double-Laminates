function [q,g,moments] = mfc_exact_spectrum(model,low,high)
%MFC_EXACT_SPECTRUM Exact paired uniform-angle moments in degrees.
aa=low*pi/180;bb=high*pi/180;width=bb-aa;
c2=(sin(2*bb)-sin(2*aa))/(2*width);c4=(sin(4*bb)-sin(4*aa))/(4*width);
m2=(1+c2)/2;n2=(1-c2)/2;
m4=(3+4*c2+c4)/8;n4=(3-4*c2+c4)/8;m2n2=(1-c4)/8;
a=model.Q(1,1);b=model.Q(2,2);c=model.Q(1,2);d=model.Q(3,3);
q11=a*m4+b*n4+2*(c+2*d)*m2n2;q22=a*n4+b*m4+2*(c+2*d)*m2n2;
q12=(a+b-4*d)*m2n2+c*(m4+n4);q66=(a+b-2*c-2*d)*m2n2+d*(m4+n4);
q=[q11 q12 0;q12 q22 0;0 0 q66];g1=model.g(1);g2=model.g(2);
g=[g1*m2+g2*n2;g1*n2+g2*m2;0];
moments=struct('cos2',c2,'cos4',c4,'variance_gx', ...
    max(0,g1*g1*m4+g2*g2*n4+2*g1*g2*m2n2-g(1)^2));
end
