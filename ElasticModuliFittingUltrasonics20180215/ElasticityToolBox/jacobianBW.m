function [J,velc]=jacobianBW(Co,iconst,sym,dcos,idfnt,rho,Cflg)
% function to calculate Jacobian of velocities wrt C (Cflg=1) or S
% J is npts by number of constants
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

npts=length(idfnt);
ncj=length(iconst);
J=zeros(npts,ncj);
deltC=1e-3*ones(ncj,1);   %increment by 0.001 GPa for finite differences
deltS=1e-6*ones(ncj,1);
So=Ci2Cij(inv(Ci2Cij(Co,sym)),sym);
velc=xstl(dcos,rho,Ci2Cij(Co,sym));
jj=1;
%Calculate matrix of partials:
for gnct=1:ncj, 
    if Cflg
      Co(iconst(gnct))=Co(iconst(gnct))+deltC(gnct); % calculate upper perturbed values 
      zH=xstl(dcos,rho,Ci2Cij(Co,sym));
      Co(iconst(gnct))=Co(iconst(gnct))-deltC(gnct); % Return to reference moduli
      J(:,jj)=(zH(idfnt)-velc(idfnt))/deltC(gnct);  % finite diference over + range 
    else
      So(iconst(gnct))=So(iconst(gnct))+deltS(gnct); % calculate upper perturbed values  
      zH=xstl(dcos,rho,inv(Ci2Cij(So,sym)));
      So(iconst(gnct))=So(iconst(gnct))-deltS(gnct); % Return to reference moduli
      J(:,jj)=(zH(idfnt)-velc(idfnt))/deltS(gnct);  % finite diference over + range 
    end    
    jj=jj+1;
end