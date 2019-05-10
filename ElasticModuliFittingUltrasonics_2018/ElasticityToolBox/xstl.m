function [velocities,eigvec] = xstl(dcos,rho,C)
% Function to calculate velocities for arbitrary symmetry and propagation direction
%  Usage:
%     ([P Sfast Sslow], eigenvectors) = xstl(dcos,rho,C)
% dcos are direction cosines
% rho is density 
% C is matrix of elastic moduli
% units: GPa for moduli, km/s for velocities, gm/cc for density
% 
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

[np,nv]=size(dcos);
if nv~=3,error('dcos needs to be number of points by 3 matrix'),end
[n,m]=size(C);
if (n~=6 ||m~=6),error('moduli need to be in 6x6 matrix'),end
if min(eig(C))<0,error('elastic moduli need to be positive definite'),end

%set up output matrixes:
velocities=zeros(np,3);
eigvec=zeros(3,3,np);
%Ct=Tnsr2Mtrx(C); % convert matix to tensor version
%loop over all directions 
A=zeros(3,3);
for iang=1:np
% Here is the way the Christoffel matrix is set up
%             A(r,s)=C(r,l,s,m)n(l)n(m) where r,s,l,m=1:3
% the derivation of this starts with the wave equation:
% rho d2ur/dt2 = Crlsm d2us/dxldxm
% substituting in u=exp(i(kx-wt))
% gives: (Crlsm kl km - rho w^2 delta(rs)) = 0
% here is the tensor form expanded out
% for r=1:3
%     for s=1:3
%         for l=1:3
%              for m=1:3
%                  A(r,s)=A(r,s)+Ct(r,l,s,m)*dcos(iang,l)*dcos(iang,m);
%              end
%         end
%     end
% end
% And below is faster code to do the same thing.
A(1,1)=C(1,1)*dcos(iang,1)^2 +C(6,6)*dcos(iang,2)^2 +C(5,5)*dcos(iang,3)^2+2*C(5,6)*dcos(iang,2)*dcos(iang,3)+2*C(1,5)*dcos(iang,3)*dcos(iang,1)+2*C(1,6)*dcos(iang,1)*dcos(iang,2);
A(2,2)=C(6,6)*dcos(iang,1)^2 +C(2,2)*dcos(iang,2)^2 +C(4,4)*dcos(iang,3)^2+2*C(2,4)*dcos(iang,2)*dcos(iang,3)+2*C(4,6)*dcos(iang,3)*dcos(iang,1)+2*C(2,6)*dcos(iang,1)*dcos(iang,2);
A(3,3)=C(5,5)*dcos(iang,1)^2 +C(4,4)*dcos(iang,2)^2 +C(3,3)*dcos(iang,3)^2+2*C(3,4)*dcos(iang,2)*dcos(iang,3)+2*C(3,5)*dcos(iang,3)*dcos(iang,1)+2*C(4,5)*dcos(iang,1)*dcos(iang,2);
A(2,3)=C(5,6)*dcos(iang,1)^2 +C(2,4)*dcos(iang,2)^2 +C(3,4)*dcos(iang,3)^2+(C(2,3)+C(4,4))*dcos(iang,2)*dcos(iang,3)+(C(3,6)+C(4,5))*dcos(iang,3)*dcos(iang,1)+(C(2,5)+C(4,6))*dcos(iang,1)*dcos(iang,2);
A(1,3)=C(1,5)*dcos(iang,1)^2 +C(4,6)*dcos(iang,2)^2 +C(3,5)*dcos(iang,3)^2+(C(3,6)+C(4,5))*dcos(iang,2)*dcos(iang,3)+(C(1,3)+C(5,5))*dcos(iang,3)*dcos(iang,1)+(C(1,4)+C(5,6))*dcos(iang,1)*dcos(iang,2);
A(1,2)=C(1,6)*dcos(iang,1)^2 +C(2,6)*dcos(iang,2)^2 +C(4,5)*dcos(iang,3)^2+(C(2,5)+C(4,6))*dcos(iang,2)*dcos(iang,3)+(C(1,4)+C(5,6))*dcos(iang,3)*dcos(iang,1)+(C(1,2)+C(6,6))*dcos(iang,1)*dcos(iang,2);
A(3,1)=A(1,3);
A(3,2)=A(2,3);
A(2,1)=A(1,2);

[v,d]=eig(A); % solve for velocities and polarizations
[vel,id]=sort(sqrt(diag(d)/rho));  %determine the order from fast to slowest velocities
velocities(iang,:)=vel(:)';  % load results in output variables
eigvec(:,:,iang)=v(:,id);
end
