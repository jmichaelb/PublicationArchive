function [chisqr,J,dv,sigdat]=jacobian_ea(InputStrct,ix,Co)
% function to calculate Jacobian of velocities for a single orientation wrt ea  
% J is npts by number of constants
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013


if(strcmp(InputStrct.Data.DataSource,'BodyWaves') || strcmp(InputStrct.Data.DataSource,'B&SWaves'))
        [veldat,sigdat,dcos,idfnt]=Data2matrixBW(InputStrct,ix);
        veldatbw=veldat(idfnt);
        sigdatbw=sigdat(idfnt);
        nptsbw=length(idfnt);
        cm=Ci2Cij(Co,InputStrct.Data.sym);
        rho=InputStrct.Data.rho;
        Jbw=zeros(nptsbw,3);
        delt_ea=1e-3*ones(3,1);   %increment by 0.001 degree for finite differences
        velR=xstl(dcos,rho,cm);
        jj=1;
        %Calculate matrix of partials:
        for gnct=1:3,
            InputStrct.opts.eulerangles(gnct,ix)=InputStrct.opts.eulerangles(gnct,ix)+delt_ea(gnct); % increment reference moduli
            [~,~,dcos,idfnt]=Data2matrixBW(InputStrct,ix);
            velH=xstl(dcos,rho,cm);
            InputStrct.opts.eulerangles(gnct,ix)=InputStrct.opts.eulerangles(gnct,ix)-delt_ea(gnct); % Return to reference moduli
            Jbw(:,jj)=(velH(idfnt)-velR(idfnt))/delt_ea(gnct);  % finite diference over + range
            jj=jj+1;
        end
        dvbw=veldatbw-velR(idfnt);
        chisqrbw=sum((dvbw./sigdatbw).^2)/nptsbw;
end
if(strcmp(InputStrct.Data.DataSource,'SurfaceWaves') || strcmp(InputStrct.Data.DataSource,'B&SWaves'))
        [veldatsw,sigdatsw,dcos]=Data2matrixSW(InputStrct,ix);
        InputStrct.opts.ifit=ix;
        nptssw=length(veldatsw);
        Jsw=zeros(nptssw,3);
        delt_ea=1e-3*ones(3,1);   %increment by 0.001 degree for finite differences
        SWvels=SurfaceWaveVel(InputStrct,Co,'v');
        velR=SWvels.velc;
        jj=1;
        %Calculate matrix of partials:
        for gnct=1:3,
            InputStrct.opts.eulerangles(gnct,ix)=InputStrct.opts.eulerangles(gnct,ix)+delt_ea(gnct); % increment reference moduli
            SWvels=SurfaceWaveVel(InputStrct,Co,'v');
            velH=SWvels.velc;
            InputStrct.opts.eulerangles(gnct,ix)=InputStrct.opts.eulerangles(gnct,ix)-delt_ea(gnct); % Return to reference moduli
            Jsw(:,jj)=(velH(:)-velR(:))/delt_ea(gnct);  % finite diference over + range
            jj=jj+1;
        end
        dvsw=veldatsw-velR;
        chisqrsw=sum((dvsw./sigdatsw).^2)/nptssw;          
end

if(strcmp(InputStrct.Data.DataSource,'B&SWaves'))
    J=[Jbw;Jsw];
    chisqr=(nptsbw*chisqrbw+nptssw*chisqrsw)/(nptsbw+nptssw);
    dv=[dvbw(:);dvsw(:)];
    sigdat=[sigdatbw(:);sigdatsw(:)];
elseif(strcmp(InputStrct.Data.DataSource,'BodyWaves'))
    J=Jbw;
    chisqr=chisqrbw;
    dv=dvbw(:);
    sigdat=sigdatbw(:);
elseif(strcmp(InputStrct.Data.DataSource,'SurfaceWaves'))
    J=Jsw;
    chisqr=chisqrsw;
    dv=dvsw(:);
    sigdat=sigdatsw(:);
end
    
