function SWPlot(Results,ifig,pltdev)
% function to plot surface wave data sets
% usage: 
%     SWPlot(Results,ifig, pltdev)
% where:
%     Results is the results structure, 
%     ifig  is the window you wish to plot in
%     pltdevs sets the percentage range for the deviation plot 
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

InputStrct=Results.InputStrct;
SWResults=Results.SWResults;
constants=InputStrct.opts.constants;
ifit=InputStrct.opts.ifit;

pltangs=0:180;
na=length(pltangs);
InputStrct.opts.ns=400;
Temp=InputStrct;
for i=ifit,
    Temp.Data.sample(i).SWangles=pltangs;
    Temp.Data.sample(i).SWivels=1:na;
end
SWout=SurfaceWaveVel(Temp,constants(:),'s');

nplt=length(ifit);
np=zeros(1,nplt);
for i=ifit
      ivels=InputStrct.Data.sample(i).SWivels;
      np(i)=length(ivels);
end
nc=[0 cumsum(np)];

if nplt==1
    ir=2;ic=1;
    iplt=1;
elseif (nplt==3)  
    ir=2;ic=3;
    iplt=[1 2 3];
elseif (nplt==4)
    ir=2;ic=4;
    iplt=[1 2 3 4];
elseif (nplt==5)
    ir=4;ic=3;
    iplt=[1 2 3 7 8];
elseif (nplt==6)
    ir=4;ic=3;
    iplt=[1 2 3 7 8 9];
elseif (nplt==7)
    ir=4;ic=4;
    iplt=[1 2 3 4 9 10 11];
elseif (nplt==8)
    ir=4;ic=4;
    iplt=[1 2 3 4 9 10 11 12];
elseif (nplt==9)
    ir=4;ic=5;
    iplt=[1 2 3 4 5 11 12 13 14];
elseif (nplt==10)
    ir=4;ic=5;
    iplt=[1 2 3 4 5 11 12 13 14 15];
end


figure(ifig)
clf
ii=0;
for i=ifit,
    ii=ii+1;
    subplot(ir,ic,iplt(ii))
    colormap('gray')
    spectra=SWout.sample(ii).spectra;
    delv=SWResults.dv(nc(i)+(1:np(i)));
    ivels=InputStrct.Data.sample(i).SWivels;
    veldata=InputStrct.Data.sample(i).SWvelocities(ivels);
    angles=InputStrct.Data.sample(i).SWangles(ivels);
    sigs=InputStrct.Data.sample(i).SWsigvels(ivels);
    txt=InputStrct.Data.sample(i).name;
    imagesc(pltangs,-spectra(:,1),log(abs(.03+spectra(:,3:end))))
    hold on
    plot(angles,-veldata,'Marker','o','Color',[0 0 0],'LineStyle','none',...
    'MarkerFaceColor',[0 0 0]) 
    %errorbar(angles,-veldata,sigs,'k.')
    ylabel('VELOCITY (km/s)')
    xlabel('ANGLE')
    title(txt)
    hold off

subplot(ir,ic,iplt(ii)+ic)
errorbar(angles,100*delv./veldata,100*sigs./veldata,'Marker','o','Color',[0 0 0],'LineStyle','none',...
    'MarkerFaceColor',[0 0 0]);
    hold on
    if(max(angles)>90),
       plot([0 180],[-.3 -.3],'k:',[0 180],[.3 .3],'k:')
       axis([ 0 180 -pltdev pltdev]);
    else
      plot([0 90],[-.02 -.02],'k:',[0 90],[.02 .02],'k:') 
      axis([ 0 90 -pltdev pltdev]);
    end
    ylabel('VELOCITY DEVIATIONS %')
    xlabel('ANGLE')
    hold off  
end


    
    
  