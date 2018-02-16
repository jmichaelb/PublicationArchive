function BWPlot(Results,ifig,pltdevs)
% function to plot bodywave data sets
% usage: 
%     BWPlot(Results,ifig, pltdev)
% where:
%     Results is the results structure, 
%     ifig  is the window you wish to plot in
%     pltdevs sets the percentage range for the deviation plot 
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013


eulerangs=Results.InputStrct.opts.eulerangles;
constants=Results.InputStrct.opts.constants;
nsamp=length(Results.InputStrct.opts.ifit);
ifit=Results.InputStrct.opts.ifit;
InputStruct=Results.InputStrct;
velplt=Results.BWResults.velplt;
velc=Results.BWResults.velc;

nplt=length(ifit);
if nplt==1
    ir=2;ic=1;
    iplt=1;
elseif (nplt==3)  
    ir=2;ic=3;
    iplt=[1 2 3];
elseif (nplt==4)
    ir=2;ic=4;
    iplt=[1 2 3 4];
end
%add more if needed

figure(ifig)
clf

npold=0;
ii=0;
temp=[];
for i=ifit,
    ii=ii+1;
subplot(ir,ic,iplt(ii))
    veldata=InputStruct.Data.sample(i).BWvelocities;
    angles=InputStruct.Data.sample(i).BWangles;
    np=length(angles);
    vc=velc((1:np)+(npold),:);
    dv=100*(veldata-vc)./veldata;
    %temp=[temp; 1e3*(veldata-vc)];
    sigs=InputStruct.Data.sample(i).BWuncertainties;
    txt=InputStruct.Data.sample(i).name;
    plot(angles,veldata(:,3),'ko',angles,veldata(:,2),'k^',angles,veldata(:,1),'k*',velplt(:,1,i),velplt(:,2:4,i),'k-')
    hold on
    ylabel('VELOCITY (km/s)')
    title(txt)
    axis([-2 182 3 9.5])
    hold off
    
subplot(ir,ic,iplt(ii)+ic);
   errorbar(angles,dv(:,1),100*sigs(:,1)./veldata(:,1),'k*')
       hold on
   errorbar(angles,dv(:,2),100*sigs(:,2)./veldata(:,2),'k^')
   errorbar(angles,dv(:,3),100*sigs(:,3)./veldata(:,3),'ko')
   plot([0 180],[-.25 -.25],'k:',[0 180],[.25 .25],'k:')
   axis([ 0 180 -pltdevs pltdevs]);
   hold off
   xlabel('ANGLE')
   ylabel('% Deviation')
   npold=np+npold;
    
end
