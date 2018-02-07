function [J, velc]=jacobianSW(Co,iconst,InputStrct,Cflg)
% function to calculate Jacobian of velocities wrt C (Cflg=1) or S
% J is npts by number of constants
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013
sym=InputStrct.Data.sym;
ncomp=InputStrct.opts.ncomp;
ncj=length(iconst);
deltC=1e-3*ones(ncj,1);   %increment by 0.001 GPa for finite differences
deltS=1e-6*ones(ncj,1);   %
SWvels=SurfaceWaveVel(InputStrct,Co,'v');
velc=SWvels.velc;
npts=length(velc);
cml=Ci2Cij(Co,sym);
sml=inv(cml);
So=Ci2Cij(sml,sym);
s1=1e3*sum(sml(1:3,:)); 
Kc=[sum(sum(cml(1:3,1:3)))/9 sum(sum(sml(1:3,1:3))).^(-1)];
J=zeros(npts+ncomp,ncj);
jj=1;
%Calculate matrix of partials:
 for gnct=1:ncj,
     if Cflg
        Co(iconst(gnct))=Co(iconst(gnct))+deltC(gnct); % calculate upper perturbed values  
        cml=Ci2Cij(Co(:),sym);
        sml=inv(cml);
        Kc2=[sum(sum(cml(1:3,1:3)))/9 sum(sum(sml(1:3,1:3))).^(-1)];
        s2=1e3*sum(sml(1:3,:));
        SWvels=SurfaceWaveVel(InputStrct,Co,'v');
        z2=SWvels.velc;
        Co(iconst(gnct))=Co(iconst(gnct))-deltC(gnct); % Return to reference moduli
        if ncomp==0,
           dc=[];
        elseif ncomp==1,
           dc=Kc2(2)-Kc(2);
        elseif ncomp==2,
       %   dc=[];
        else
          dc=(s2-s1);
        end 
        J(:,jj)=[(z2(:)-velc(:)); dc(:)]/deltC(gnct);
     else
        So(iconst(gnct))=So(iconst(gnct))+deltS(gnct); % calculate upper perturbed values  
        sml=Ci2Cij(So(:),sym);
        cml=inv(sml);
        Kc2=[sum(sum(cml(1:3,1:3)))/9 sum(sum(sml(1:3,1:3))).^(-1)];
        s2=1e3*sum(sml(1:3,:));
        SWvels=SurfaceWaveVel(InputStrct,Ci2Cij(cml,sym),'v');
        z2=SWvels.velc;
        So(iconst(gnct))=So(iconst(gnct))-deltS(gnct); % Return to reference moduli
        if ncomp==0,
           dc=[];
        elseif ncomp==1,
           dc=Kc2(2)-Kc(2);
        elseif ncomp==2,
       %   dc=[];
        else
          dc=(s2-s1);
        end 
        J(:,jj)=[(z2(:)-velc(:)); dc(:)]/deltS(gnct);
     end            
    jj=jj+1;
 end
