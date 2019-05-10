function SWout=SurfaceWaveVel(InptStrct,constants,SWflg)
%  velocities=SurfaceWaveVels((InptStrct,constants,SWflg));
% where InptStrct is the standard input structure
%       constants is the vector of moduli
%       atr are orientation matrixes for upper and lower materials
%       SWflg is 'v' to return velocities at specified directions and 's'
%       to return spectra over the entire range of velocities
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             2006

npar=[1 3];  % Green's Function tensor component to use  (1 3] or [3 3];
damping=0.01; % artificial damping to avoid singularities in the Rayleigh wave calculation.
vrange=InptStrct.opts.vrange;
cml=Ci2Cij(constants,InptStrct.Data.sym);
switch SWflg(1)
    case 'v'   % calculate velocities for the experimental angles for a particlar crystal orientation
        dvcat=[];
        sig=[];
        velc=[];
        vexpcat=[];
        if min(eig(cml))>0, 
            j=0;     
        for i=InptStrct.opts.ifit,
            ivels=InptStrct.Data.sample(i).SWivels;
            j=j+1; 
            na=length(ivels);
            OM=eiler(InptStrct.opts.eulerangles(:,i)');
            v=modevel(vrange,na,InptStrct.Data.sample(i).SWvelocities(ivels),InptStrct.Data.sample(i).SWangles(ivels),InptStrct.Data.sample(i).film.orientation, ...
                 OM,InptStrct.Data.sample(i).film.constants,cml,InptStrct.Data.sample(i).film.density, InptStrct.Data.rho,...
                  InptStrct.Data.sample(i).wavelength,InptStrct.Data.sample(i).film.thickness,InptStrct.Data.sample(i).mflg,npar);                
            SWout.sample(j).velocities=v(:); 
            vexp=InptStrct.Data.sample(i).SWvelocities(ivels);
            dv=vexp-v(:);
            id=isfinite(dv);
            dvcat=[dvcat;dv(id)];
            sig=[sig; InptStrct.Data.sample(i).SWsigvels(ivels(id))];
            velc=[velc; v(:)];
            SWout.sample(i).velc=v(:);
        end
        else
            disp('non positive definite elastic tensor')
        end
        np=length(dvcat);
        SWout.rms=sqrt(sum(dvcat.^2)/np);
        SWout.dv=dvcat;
        SWout.sig=sig;
        SWout.velc=velc;
        SWout.np=np;
        
    case 's'   % calculate spectra (surface wave intensities) as a function of direction and velocity
        j=0;
        for i=InptStrct.opts.ifit,
             ivels=InptStrct.Data.sample(i).SWivels;
             j=j+1; 
             na=length(ivels);  
             OM=eiler(InptStrct.opts.eulerangles(:,i)');
            [v,y]=modeconv(InptStrct.opts.ns,InptStrct.opts.vmin,InptStrct.opts.vmax,na, ...
                  InptStrct.Data.sample(i).SWangles(ivels),InptStrct.Data.sample(i).film.orientation, ...
                  OM,InptStrct.Data.sample(i).film.constants,cml,InptStrct.Data.sample(i).film.density,...
                  InptStrct.Data.rho,InptStrct.Data.sample(i).wavelength, ...
                  InptStrct.Data.sample(i).film.thickness,damping,InptStrct.Data.sample(i).mflg,npar);                              
              SWout.sample(j).vels=v(:);                      
              SWout.sample(j).spectra=y;
        end
end
