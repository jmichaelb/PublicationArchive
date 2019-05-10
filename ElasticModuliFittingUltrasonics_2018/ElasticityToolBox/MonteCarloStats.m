function [uncerts,Cfs,rms]=MonteCarloStats(Input,npts,Cm,fac,dispflg)
% This function calculate uncertainites based on a Monte Carlo approach.
% The idea is to use the actual distribution of propagation directions and
% polarizations and use the typical experimental uncertainties. Then calculate
% a set of synthetic velocities with similar misfit statistics using a model set of
% elastic moduli.  Repeat many times and examine the resulting spread in
% moduli.  Usage:
%         [uncerts,Cfs,rms]=MonteCarloStats(Input,npts,Cm,fac,dispflg)
% where: Input is the standard structure with data for a sample
%        npts is the number of sythetic data sets to try 
%        Cm is the model set of moduli
%        fac is a multiplication factor to change the variance of the
%        synthetic data. fac=1 uses the experimental uncertainties.
%        dispflg=1 displays the fitting iterations for each model data set
%        uncert is the list of 2*the standard deviations of each modulus
%        Cfs is a nptsx#moduli matrix of all model moduli
%        rms is a vector of all model misfits.
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

nconst=length(Cm);
rms=zeros(npts,1);
Cfs=zeros(npts,nconst);
uncerts=zeros(nconst,1);
ifit=Input.opts.ifit;
if(isfield(Input.opts,'eulerangles'))
 ea=Input.opts.eulerangles;
 Input.Data.dcosflg=0;
else
    Input.Data.dcosflg=1;
    ea=[];
end
sym=Input.Data.sym;
rho=Input.Data.rho;

if (isfield(Input.Data.sample,'SWvelocities'))
    if not(isfield(Input.Data.sample,'BWvelocities'))
        DataSource='SurfaceWaves';
    else
        DataSource='B&SWaves';
    end
else
    DataSource='BodyWaves';
end

for j=1:npts
switch DataSource
    case 'BodyWaves' 
        for i=ifit
            [veldat,sigdat,dcos,idfnt]=Data2matrixBW(Input,i);
            [np,~]=size(veldat);
            BWvelocities=(1+(fac*sigdat./veldat).*randn(np,3)).*xstl(dcos,rho,Ci2Cij(Cm,sym));
            velocities(idfnt)=nan*ones(size(idfnt));
            Input.Data.sample(i).BWvelocities=BWvelocities;
        end
    case 'SurfaceWaves'
         for i=ifit
            [veldat,sigdat,dcos]=Data2matrixSW(Input,i);
            np=length(veldat);
            Input.opts.ifit=i;
            SWvels=SurfaceWaveVel(Input,Cm,'v') ;
            SWvelocities=(1+(fac*sigdat./veldat).*randn(np,1)).*SWvels.velc;
            Input.Data.sample(i).SWvelocities=SWvelocities;
         end
    case 'B&SWaves'
        for i=ifit
            [veldat,sigdat,dcos,idfnt]=Data2matrixBW(Input,i);
            [np,~]=size(veldat);
            BWvelocities=(1+(fac*sigdat./veldat).*randn(np,3)).*xstl(dcos,rho,Ci2Cij(Cm,sym));
            velocities(idfnt)=nan*ones(size(idfnt));
            Input.Data.sample(i).BWvelocities=BWvelocities;
            [veldat,sigdat,dcos]=Data2matrixSW(Input,i);
            np=length(veldat);
            Input.opts.ifit=i;
            SWvels=SurfaceWaveVel(Input,Cm,'v') ;
            SWvelocities=(1+(fac*sigdat./veldat).*randn(np,1)).*SWvels.velc;
            Input.Data.sample(i).SWvelocities=SWvelocities;
        end

end
Input.opts.ifit=ifit;
Cs=Cm+3*randn(length(Cm),1);
[Cf,~,ResultsStrct]=Velocities2Cij(Input,Cs,'y',ea,'n','LM',dispflg);
rms(j)=ResultsStrct.rms;
Cfs(j,:)=Cf(:)';
end

for i=1:nconst
uncerts(i)=2*std(Cfs(:,i));
end

uncerts=round(10*uncerts)/10;