function [veldat,sigdat,dcos,idfnt]=Data2matrixBW(InptStrct,ifit)
 % the data structure is organized into individual samples.  Fitting
 % requires all data to be organized into a single list of directions and
 % velocity measurements.  If input is in terms of lab angles associated
 % with an orientation of each sample, then direction cosines need to be
 % calculated
 %
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

nfit=length(ifit);
nd=zeros(nfit,1);
for i=1:nfit
    nd(i)=length(InptStrct.Data.sample(ifit(i)).BWvelocities(:,1));
end
nsbw=[0 ; cumsum(nd)];
veldat=zeros(nsbw(end),3);
sigdat=veldat;
dcos=veldat;
for i=1:nfit,
    a_bw=InptStrct.Data.sample(ifit(i)).BWangles;
    veldat(nsbw(i)+1:nsbw(i+1),:)=InptStrct.Data.sample(ifit(i)).BWvelocities;
    sigdat(nsbw(i)+1:nsbw(i+1),:)=InptStrct.Data.sample(ifit(i)).BWuncertainties;
    if InptStrct.Data.dcosflg==1
       dcos(nsbw(i)+1:nsbw(i+1),:)=InptStrct.Data.sample(ifit(i)).dcos;
    else  
       dcos(nsbw(i)+1:nsbw(i+1),:)=angles2dcos(a_bw,InptStrct.opts.eulerangles(:,ifit(i)));
    end
end  
idfnt=find(isfinite(veldat(:)));
