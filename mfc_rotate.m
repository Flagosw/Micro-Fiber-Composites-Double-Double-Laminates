function [qb,gb] = mfc_rotate(q,g,theta)
%MFC_ROTATE Rotate orthotropic reduced tensors using engineering shear.
%   Output matrix/vector components occupy dimensions 3/4 and 3,
%   respectively. Dimensions 1/2 preserve the supplied angular grid.
%   Use cos(theta*pi/180), matching the Python numerical convention.
m=cos(theta*pi/180);n=sin(theta*pi/180);
m2=m.^2;n2=n.^2;m4=m2.^2;n4=n2.^2;mn2=m2.*n2;
a=q(1,1);b=q(2,2);c=q(1,2);d=q(3,3);
qb=zeros([size(theta),3,3]);gb=zeros([size(theta),3]);
qb(:,:,1,1)=a*m4+b*n4+2*(c+2*d)*mn2;
qb(:,:,2,2)=a*n4+b*m4+2*(c+2*d)*mn2;
qb(:,:,1,2)=(a+b-4*d)*mn2+c*(m4+n4);qb(:,:,2,1)=qb(:,:,1,2);
qb(:,:,1,3)=(a-c-2*d)*m2.*m.*n-(b-c-2*d)*m.*n2.*n;qb(:,:,3,1)=qb(:,:,1,3);
qb(:,:,2,3)=(a-c-2*d)*m.*n2.*n-(b-c-2*d)*m2.*m.*n;qb(:,:,3,2)=qb(:,:,2,3);
qb(:,:,3,3)=(a+b-2*c-2*d)*mn2+d*(m4+n4);
gb(:,:,1)=g(1)*m2+g(2)*n2;gb(:,:,2)=g(1)*n2+g(2)*m2;
gb(:,:,3)=(g(1)-g(2))*m.*n;
end
