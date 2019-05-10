function [veldat,sigdat,dcos,comp,dcomp]=Data2matrixSW(InptStrct,ifit)
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
    nd(i)=length(InptStrct.Data.sample(ifit(i)).SWangles(:));
end
nssw=[0 ; cumsum(nd)];
veldat=zeros(nssw(end),1);
sigdat=veldat;
dcos=zeros(nssw(end),3);
for i=1:nfit,
    a_sw=InptStrct.Data.sample(ifit(i)).SWangles;
    veldat(nssw(i)+1:nssw(i+1),:)=InptStrct.Data.sample(ifit(i)).SWvelocities;
    sigdat(nssw(i)+1:nssw(i+1),:)=InptStrct.Data.sample(ifit(i)).SWsigvels;
    if InptStrct.Data.dcosflg==1
       dcos(nssw(i)+1:nssw(i+1),:)=InptStrct.Data.sample(ifit(i)).dcos;
    else    
       dcos(nssw(i)+1:nssw(i+1),:)=angles2dcos(a_sw,InptStrct.opts.eulerangles(:,InptStrct.opts.ifit(i)));
    end
end  
if (isfield(InptStrct.Data,'compliances'))
  temp=InptStrct.Data.compliances;
end
ncomp=0;
comp=[];
dcomp=[];
if InptStrct.opts.constrflg(1)=='y',
   ncomp=length(temp)/2;
   if ncomp==1,
      comp=temp(1);
      dcomp=temp(2);
   elseif ncomp==3,
      comp=temp(1:3);
      dcomp=temp(4:6);
   elseif ncomp==6,
      comp=temp(1:6);
      dcomp=temp(7:12);
   else
      error('wrong number of compliance constraints')
   end
end 


