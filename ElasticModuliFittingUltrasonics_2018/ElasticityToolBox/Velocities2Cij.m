function [Cout,eaout,ResultsStrct,Ctemp]=...
               Velocities2Cij(InputStrct,Cin,Refine_Cij,ea,Refine_ea,Method,iter_disp)
%function to optimize elastic moduli based on measured velocities 
% Usage: 
%         [Cout,eaout,ResultsStrct]=Velocities2Cij(InputStrct,Cin,Refine_Cij,ea,Refine_ea,Method)
%
% where:   Cout is vector of improved moduli
%          ResultsStrct is a structure containing all the information about the fit
%          InputStrct is structure with data and fitting options
%          Cin is vector of initial guesses for elastic moduli
%          Refine_Cij is 'y' 'n' or 'r' to initialize with random values
%          ea is the matrix of euler angles (column of 3 angles for each cut)
%          Refine_ea  set to 'y' to start with current ea, set to 'n' to
%              fix current ea, and 'r' to optimize from random starting point
%          Method is  'LM' or 'NM' or 'BG'
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

DataSource='SurfaceWaves'; % Need to initialize
if nargin==6
   iter_disp=1;
   checkinput;  %(InputStrct,Cin,Refine_Cij,ea,Refine_ea,Method);
elseif nargin==7
   checkinput;  %(InputStrct,Cin,Refine_Cij,ea,Refine_ea,Method);
else
    error('incorrect number of input variables')
end

%***********************************************************************************
% use the .opts field of the data structure to hold (temporary and output) information 
% about the current optimization attempt.  These values will change 
%***********************************************************************************
InputStrct.opts.eulerangles=ea;  
InputStrct.opts.constants=Cin;
iconst=InputStrct.opts.iconst;  % if some constants are to be left out of the fitting
nconst=length(iconst);
ifit=InputStrct.opts.ifit;
sym=InputStrct.Data.sym;
rho=InputStrct.Data.rho;
% set a few other things
lb=InputStrct.Data.Trust.constants(:,1);
ub=InputStrct.Data.Trust.constants(:,2);
ub=ub(iconst);
lb=lb(iconst);
Cout=Cin;   % default is to return the input values
eaout=ea;
Ctemp=Cout;
time=0;
switch Refine_Cij(1) 
    case 'y'
        Co=Cin(iconst);  % Co has just the moduli being adjusted
    case 'r'  %start with random moduli
        Ctemp=Crand(InputStrct.Data.Trust.constants); 
        Co=Ctemp(iconst); % set the moduli to vary to random values (within trust region)
        Ctemp=Co;  % this sets the other moduli equal to the input
    case 'n'
        Co=Cin;   %  Co is full length
end

%***********************************************************************************
% determine which calculations are needed and organize data:
%***********************************************************************************
switch DataSource
    case 'BodyWaves'
        rmsfunc='BW_calc';
        [veldat_bw,sigdat_bw,dcos_bw,idfnt]=Data2matrixBW(InputStrct,ifit);
        npts_bw=length(idfnt);             
    case 'SurfaceWaves'
        rmsfunc='SW_calc';
        [veldat_sw,sigdat_sw,~,comp,dcomp]=Data2matrixSW(InputStrct,ifit);
        ncomp=length(dcomp);
        InputStrct.opts.ncomp=ncomp;
        npts_sw=length(veldat_sw);
    case 'B&SWaves'
        rmsfunc='EC_calc';
        [veldat_bw,sigdat_bw,dcos_bw,idfnt]=Data2matrixBW(InputStrct,ifit);
        npts_bw=length(idfnt);
        [veldat_sw,sigdat_sw,~,comp,dcomp]=Data2matrixSW(InputStrct,ifit);
        ncomp=length(dcomp);
        npts_sw=length(veldat_sw);
        InputStrct.opts.ncomp=ncomp;
end

%***********************************************************************************
% set up the fitting method
%***********************************************************************************
switch Method
    case 'NM' % Nelder-Mead Method
        options=optimset('Display','off','MaxIter',700);    
        functxt=['[Co,newmisfit,~,output]=fminsearch(@' rmsfunc ',Co,options);'];
    case 'BG'  % Backus-Gilbert Inverse
        functxt='[Co,newmisfit]=BackusGilbert(Co);';
    case 'LM' % Levenberg-Marquardt Inverse
        functxt='[Co,newmisfit]=LM_LSQR(Co);';  
        InputStrct.opts.funiter=0;
end

%********************************************************************************************
% Optimize the Cij
%***********************************************************************************
if (strcmp(Refine_Cij(1),'y') || strcmp(Refine_Cij(1),'r') && strcmp(Refine_ea,'n')) 
   if(iter_disp),tic,end 
   oldmisfit=10000;
   if (strcmp(Method,'LM') )
       txt='iteration    chisqr      optimality       lambda      relaxation';
       if(iter_disp),disp(txt),end
   elseif (strcmp(Method,'BG'))
       txt='iteration    chisqr      optimality      variance     relaxation';
       if(iter_disp),disp(txt),end
   elseif (strcmp(Method,'NM'))
       txt='iteration    chisqr  ';
       if(iter_disp),disp(txt),end
       initial=0;
       eval(['initial=' rmsfunc '(Co);'])
       txt=sprintf('%4i        %8.2f   ',0,initial);
       if(iter_disp),disp(txt),end
   end 
   iter=0;
   itrial=0;
   while (oldmisfit>150 && iter<100 && itrial<3)  % outer loop to try again if fitter stalls
     testmf=1; 
     count=0;
     itrial=itrial+1;
     output=[];
     newmisfit=0; 
   while (testmf>1e-3 && count<10)  % inner loop to run optimizer
        count=count+1;
        eval(functxt)
        if strcmp(Method,'NM')
         iter=iter+output.iterations;
         txt=sprintf('%4i        %8.2f   ',iter,newmisfit);
         if(iter_disp),disp(txt),end
        end
        testmf=abs(oldmisfit-newmisfit)/oldmisfit;
        oldmisfit=newmisfit;
   end
   end
Cout(iconst)=Co;
InputStrct.opts.constants=Cout;
if(iter_disp),time=toc;else,time=0;end

%***********************************************************************************
%Optimize euler angles
%***********************************************************************************
elseif(strcmp(Refine_Cij(1),'n') && (strcmp(Refine_ea,'y') || strcmp(Refine_ea,'r'))) 
  tic
  lbea=ea-InputStrct.Data.Trust.eulerangles;
  ubea=ea+InputStrct.Data.Trust.eulerangles;
  for ie=ifit       
     ea_temp=ea(:,ie);
     [eaout(:,ie),chisqr]=LM_LSQR_ea(ea_temp,ie,lbea(:,ie),ubea(:,ie));
     if strcmp(Refine_ea(1),'r')
            icount=0;
            best=chisqr;
            testea=-100;
            countmax=100;
            while (testea<1e-6 && icount<countmax)   
                 ea_temp=lbea(:,ie)+(ubea(:,ie)-lbea(:,ie)).*rand(3,1);
                [eatrial,new]=LM_LSQR_ea(ea_temp,ie,lbea(:,ie),ubea(:,ie));
                testea=(best-new)/best;
                icount=icount+1;
            end
            if icount<countmax  % only change euler angles if loop finished early with a better solution
                eaout(:,ie)=eatrial; %
            end
     end
  end
InputStrct.opts.eulerangles=eaout;
eachange=max(max(abs(ea-eaout)./InputStrct.Data.Trust.eulerangles));
if eachange>0.9
        disp('Adjusted eulerangles may be close to trust boundary. Is an increase in the trust region needed?')
end
if strcmp(DataSource, 'BodyWaves')
    [~,~,dcos_bw]=Data2matrixBW(InputStrct,ifit);
elseif strcmp(DataSource, 'SurfaceWaves')
    [~,sigdat_sw,~,comp,dcomp]=Data2matrixSW(InputStrct,ifit);
end
ea=eaout; %reset the euler angles for further calculations
time=toc;
end

%***************************************************************************************************
% Determine statistics related to current full set of moduli and put
% results in output structure
%***********************************************************************************
nconst=length(Cout);
iconst=1:nconst;
cmo=Ci2Cij(Cout,sym);
smo=inv(cmo);
So=Ci2Cij(smo,sym);
s=1e3*sum(smo(1:3,:));
ResultsStrct.datestamp=datestr(now);
ResultsStrct.InputStrct=InputStrct;

if (strcmp(DataSource,'BodyWaves') || strcmp(DataSource, 'B&SWaves'))
    [JCbw,~]=jacobianBW(Cout,1:nconst,sym,dcos_bw,idfnt,rho,1);
    [JSbw,velcbw]=jacobianBW(Cout,1:nconst,sym,dcos_bw,idfnt,rho,0); 
    [chisqrbw,~,dv_bw,rmsbw,~]=BW_calc(Cout);
    sigexbw=abs(dv_bw);
    ids=find(sigexbw-sigdat_bw(idfnt)<0);
    sigexbw(ids)=sigdat_bw(idfnt(ids)); % individual data uncertainties are no smaller than assumed uncertainty
    rms=rmsbw;
    chisq=chisqrbw;
    Acbw=JCbw./(sigexbw*ones(1,nconst));
    Bcbw=inv(Acbw'*Acbw);
    Two_sigCbw=2*sqrt(diag(Bcbw)); 
    Asbw=JSbw./(sigexbw*ones(1,nconst));
    Bsbw=inv(Asbw'*Asbw);
    Two_sigSbw=1e3*2*sqrt(diag(Bsbw));
    out=KG_calc(Cout,Bcbw,Bsbw,sym);
    ResultsStrct.CijSij=round(10*[Cout(:) Two_sigCbw(:) 1e3*So(:) Two_sigSbw(:)])/10;
            
%calculate points for body wave plots  
if not(isempty(ea))
    angles_plt=0:InputStrct.opts.pltRange;
    na=length(angles_plt);
    dcos_plt=zeros(na,3,length(ifit));
    BWvel_plt=zeros(na,4,length(ifit));
    for iplt=ifit
        dcos_plt(:,:,iplt)=angles2dcos(angles_plt,ea(:,iplt));
        BWvel_plt(:,2:4,iplt)=xstl(squeeze(dcos_plt(:,:,iplt)),rho,Ci2Cij(Cout,sym));
        BWvel_plt(:,1,iplt)=angles_plt(:);
    end 
    BWResults.velplt=BWvel_plt;
end    
     BWResults.Two_sig.C=Two_sigCbw;
     BWResults.Two_sig.S=Two_sigSbw;
     BWResults.constants=Cout;
     BWResults.rms=rmsbw;
     BWResults.chisq=chisq;
     BWResults.velc=velcbw;
     BWResults.sigdat=sigdat_bw;
     BWResults.dv=dv_bw; 
     BWResults.idfnt=idfnt; 
     BWResults.JC=JCbw;
     BWResults.JS=JSbw;
     BWResults.covariance=Bcbw;
     ResultsStrct.BWResults=BWResults; 
     ResultsStrct.CijSij=[round(10*Cout(:))/10 round(10*Two_sigCbw(:))/10 ...
        round(1e5*[So(:) 1e-3*Two_sigSbw(:)])/100];

end

if (strcmp(DataSource,'SurfaceWaves') || strcmp(DataSource, 'B&SWaves'))
    [JCsw, ~]=jacobianSW(Cout,1:nconst,InputStrct,1);
    [JSsw, ~]=jacobianSW(Cout,1:nconst,InputStrct,0);
    [chisqrsw,~,dv_sw,rmssw,~]=SW_calc(Cout);
    SWvels=SurfaceWaveVel(InputStrct,Cout,'v') ;
    Kc=[sum(sum(cmo(1:3,1:3)))/9 sum(sum(smo(1:3,1:3))).^(-1)];
    sigexsw=abs(dv_sw);
    rms=rmssw;
    ids=find((sigexsw(1:npts_sw)-sigdat_sw)<0);
    sigexsw(ids)=sigdat_sw(ids);
    if ncomp==0
       dc=[];
    elseif ncomp==1
      dc=comp-Kc(2);
    elseif ncomp==2
      dc=comp-a;   % not implemented correctly
    else
      dc=(comp-s);
    end 
    SWResults.chicomp=dc;  
    idc=find((dcomp-abs(dc))<0);
    temp=dcomp;
    temp(idc)=dc(idc);
    sigexsw=[sigexsw(1:npts_sw)' temp(:)'];
    chisq=chisqrsw;
    Acsw=JCsw./(sigexsw(:)*ones(1,nconst));
    Bcsw=inv(Acsw'*Acsw);
    Two_sigCsw=2*sqrt(diag(Bcsw));
    Assw=JSsw./(sigexsw(:)*ones(1,nconst));
    Bssw=inv(Assw'*Assw);
    Two_sigSsw=1e3*2*sqrt(diag(Bssw));
    SWResults.Two_sig.C=Two_sigCsw;
    SWResults.Two_sig.S=Two_sigSsw;
    SWResults.rms=rmssw;
    SWResults.dv=SWvels.dv; 
    ResultsStrct.SWResults=SWResults;
    ResultsStrct.SWResults.sigex=sigexsw(:);
    out=KG_calc(Cout,Bcsw,Bssw,sym);
    ResultsStrct.CijSij=[round(10*Cout(:))/10 round(10*Two_sigCsw(:))/10 ...
        round(1e5*[So(:) Two_sigSsw(:)/1e3])/100];
end  

if strcmp(DataSource,'B&SWaves')
    AJC=[JCbw;JCsw]./([sigexbw(:);sigexsw(:)]*ones(1,nconst));        
    BJc=inv(AJC'*AJC);
    Two_sigCJ=2*sqrt(diag(BJc));
    rms=sqrt((npts_bw*rmsbw^2+npts_sw*rmssw^2)/(npts_bw+npts_sw));
    chisq=(npts_bw*chisqrbw+(npts_sw+ncomp)*chisqrsw)/(npts_bw+npts_sw+ncomp);
    AJS=[JSbw;JSsw]./([sigexbw(:);sigexsw(:)]*ones(1,nconst));        
    BJs=inv(AJS'*AJS);
    Two_sigSJ=2*sqrt(diag(BJs)) ;
    out=KG_calc(Cout,BJc,BJs,sym);
    ResultsStrct.CijSij=[round(10*Cout(:))/10 round(10*Two_sigCJ(:))/10 ...
        round(1e5*[So(:) Two_sigSJ(:)])/100];
end

sij=Ci2Cij(So,sym);

ResultsStrct.rms=rms*1e3;
ResultsStrct.chisq=chisq;
ResultsStrct.KG=out;
ResultsStrct.comp=round(sum(sij(1:3,:))*1e5)/1e2;
ResultsStrct.Method=Method;
ResultsStrct.RefineCij=Refine_Cij;
ResultsStrct.Refineea=Refine_ea;

txt=sprintf('rms misfit =%4.1f m/s  chisqr = %6.2f  elapsed time %4.1f s \n',...
    rms*1e3,chisq, time);
if(iter_disp),disp(txt),end
%*************************************************************************************************
% end of main routine
%*************************************************************************************************


%***********************************************************************************
% Below here are nested functions used in the optimization - these share  variables (shown
% in cyan) with main function
%***********************************************************************************

function  [chisqr,J,dvbw,rms,npflg] = BW_calc(parms)
% Function that calculates the misfit of bodywave data for a given set of moduli
% Usage:
%    [chsqr,J,dvbw,rms,npflg] = BW_calc(parms)
%  where:
%            parms is the vector of elastic moduli 
%            chisqr is the sum of square of deviations between observations and the 
%                 predictions weighted by uncertainties
%            rms is the misfit without uncertainty weighting
%            J is the Jacobian
%            dvbs is vector of velocity deviations
%            npflg is set to 1 if the moduli are not positive definite
%  calls:    eiler (to turn the euler angles into an orientation matrix)
%            angles2dcos to turn lab angles into direction cosines
%            Ci2Cij (to turn the vector constants into a 6x6 matrix)
%            xstl (to calculate velocities for the given directions and constants)
%            jacobianBW to calculate partial derivatives of velocities wrt moduli
%   needs:   InputStrct structure with all the data
% 
%  jmb 5/2010  8/2013

% work with local variable for moduli and direction cosines - don't change the global values
C_bw=InputStrct.opts.constants; 
dcos_tbw=dcos_bw;   
C_bw(iconst)=parms;            % load moduli that are being adjusted
c=Ci2Cij(C_bw,sym);
if (min(eig(c))>0)
      npflg=0;
      velcc=xstl(dcos_tbw,rho,c);
      dvbw=veldat_bw(idfnt)-velcc(idfnt);
      sigexp=sigdat_bw(idfnt); 
      rms=sqrt(sum(dvbw.^2)/npts_bw);
      chisqr=sum((dvbw./sigexp).^2)/npts_bw;  
      if not(strcmp(Method,'NM'))          
         J=jacobianBW(C_bw,iconst,sym,dcos_tbw,idfnt,rho,1);
      else
          J=[];
      end
else
      npflg=1;
      np=length(idfnt);
      dvbw=ones(np,1);
      rms=1e6;
      chisqr=1e6;
      J=ones(np,length(iconst)); 
end  %check of positive definite  
end

%***********************************************************************************
function  [misfit,J,dv,rms,npflg] = SW_calc(parms)
% Function that calculates the misfit of surfacewave data for a given set of moduli
% Usage:
%    [chsqr,J,dvbw,rms,npflg] = SW_calc(parms)
%  where:
%            parms is the vector of elastic moduli  
%            chisqr is the sum of square of deviations between observations and the 
%                 predictions weighted by uncertainties
%            rms is the velocity misfit without uncertainty weighting
%            J is the Jacobian
%            dvbs is vector of velocity deviations
%            npflg is set to 1 if the moduli are not positive definite
%  calls:    SurfaceWaveVel to calculate surface wave velocities
%            Ci2Cij (to turn the vector constants into a 6x6 matrix)
%            jacobianSW to calculate partial derivatives of velocities wrt
%            moduli
%   needs:   InputStrct structure with all the data
% 
%  jmb 5/2010  8/2013
C_sw=InputStrct.opts.constants;
C_sw(iconst)=parms;            % load moduli that are being adjusted
cm=Ci2Cij(C_sw,InputStrct.Data.sym);
sm=inv(cm);
sl=1e3*sum(sm(1:3,:));
K=[sum(sum(cm(1:3,1:3)))/9 sum(sum(sm(1:3,1:3))).^(-1)];  % Voigt-Reuss Bulk Modulus
if InputStrct.opts.constrflg(1)=='y'
        if ncomp==0
            delconstr=[];
        elseif ncomp==1
            delconstr=(comp-K(2))/dcomp;
        elseif ncomp==2
            % dc=[];
        elseif ncomp==6
            delconstr=(comp-sl)./dcomp;      
        else
            error('wrong number of constraints')
        end 
else
        delconstr=[];
end

if (min(min(isfinite(cm)))) 
    if (min(eig(cm))>0)
        npflg=0;
        SWvels=SurfaceWaveVel(InputStrct,C_sw,'v')  ;   
        rms=SWvels.rms;  
        dvsw=SWvels.dv;
        misfit=[(dvsw./sigdat_sw); delconstr(:)];
        dv=[dvsw;delconstr(:).*dcomp(:)];
        chisq=sum(misfit.^2)/(npts_sw+ncomp) ;
        switch Method
            case 'NM'
                misfit=chisq;
                J=[];
            case {'LM' , 'BG'}
                J=jacobianSW(C_sw,iconst,InputStrct,1);
                misfit=chisq;
        end
    else
        npflg=1;
        nd=length(sigdat_sw);
        misfit=1e6;
        J=ones(nd+ncomp,nconst);
        dv=ones(nd+ncomp,1);
        rms=1e6;
    end
else    % failure - set large misfit
        npflg=1;
        nd=length(sigdat_sw);
        misfit=1e6;
        J=ones(nd+ncomp,nconst);
        dv=ones(nd+ncomp,1);
        rms=1e6;
end
end
%***********************************************************************************
function [chisqr,J,dv,rms,npflg] = EC_calc(parms)
% Function that calculates the misfit of surfacewave and bodywave data for a given set 
% of moduli
% Usage:
%    [chsqr,J,dvbw,rms,npflg] = EC_calc(parms)
% Calls BW_calc and SW_calc

    [chisqrBW,JBW,dvBW,rmsBW,npflg] = BW_calc(parms);
    [chisqrSW,JSW,dvSW,rmsSW,~] = SW_calc(parms);
    chisqr=((npts_sw+ncomp)*chisqrSW  + npts_bw*chisqrBW)/(npts_bw+npts_sw+ncomp);
    rms=  sqrt((npts_sw*rmsSW^2  + npts_bw*rmsBW^2)/(npts_bw+npts_sw));
    dv=[dvBW;dvSW];
switch Method
    case 'NM'
         J=[];
    case {'LM' , 'BG'}
         J=[JBW;JSW];
end
end

%***********************************************************************************
function [Cout,chisqr]=BackusGilbert(Cin)
%function to optimize elastic moduli using the inverse method of Backus and Gilbert as
%   described in Weidner & Carleton 1977
% Usage:   Cout=BackusGilbert(Cin)
% where Cin and Cout are vectors of elastic constants
% InputStrct must be an avaliable global structure for this routine to work
% calls BW_calc, SW_calc, or EC_calc
%jmb 3/2010
testcw=1000;
relax=.3;
oldchisqr=1e6;
nc=length(Cin);
Cout=Cin;
i_iter=InputStrct.opts.funiter;
i_stop=i_iter+10;
%flglam=1;
istuck=0;
while (abs(testcw)>1e-4 && i_iter<i_stop && relax>1e-3 && istuck<4)
    switch DataSource
        case 'BodyWaves'
            [chisqr,J,dv,~,npflg] = BW_calc(Cout);
            J=J';
        case 'SurfaceWaves'
            [chisqr,J,dv,~,npflg]= SW_calc(Cout);
            J=J';
        case 'B&SWaves'
           [chisqr,J,dv,~,npflg] = EC_calc(Cout);
            J=J';
    end
    
   ndat=length(dv);
   vr=var(dv);
   alpha=zeros(nc,ndat);
   for i=1:nc  % solve for alpha row by row
     idd=find((1:13)~=i);  
     Amj=J(idd,1:ndat)'*J(idd,1:ndat) + vr*eye(ndat); 
     A=[[Amj J(i,1:ndat)'];[J(i,1:ndat) 0]];
     d=[zeros(ndat,1); 1];
     sol=A\d;
     alpha(i,:)=sol(1:ndat)';
   end  
   dC=alpha*(dv(:));
   testcw=(oldchisqr-chisqr)/chisqr;
   if (testcw>0 && npflg==0)
        txt=sprintf('%4i       %8.3f      %8.3e       %8.3e     %8.3e  ',...
            i_iter+InputStrct.opts.funiter,chisqr,testcw,vr,relax);
        if(iter_disp),disp(txt),end
        i_iter=i_iter+1;
        Cbest=Cout;
        Cout=Cout(:)+relax*dC;
        oldchisqr=chisqr;
        %lam=lam/5;
        relax=relax*1.25;
        if relax>4,relax=4;end
    elseif (testcw<=0 && npflg==0)
        relax=relax/1.25;
        Cout=Cout(:)+relax*dC; 
        istuck=istuck+1;
    elseif (npflg==1)
        %lam=lam*10; 
        relax=relax/10; 
        istuck=istuck+1;
    end
    id=find((Cout-lb)<0);
    Cout(id)=lb(id);
    id=find((Cout-ub)>0);
    Cout(id)=ub(id); 
end
Cout=Cbest;
InputStrct.opts.funiter=i_iter;
end

%***********************************************************************************
function [Cout,chisqr]=LM_LSQR(Cin)
%function to optimix elastic constants using the Levenberg-Marquardt method
% Usage:   Cout=LM_LSQR(Cin)
% where Cin and Cout are vectors of elastic constants
% InputStrct must be an avaliable global structure for this routine to work
% calls BW_calc
%jmb 3/2010
relax=1e0;
Cout=Cin;
Cbest=Cin;
testcw=1000;
lam=1e-2;
oldchisqr=1e6;
i_iter=InputStrct.opts.funiter;
i_stop=i_iter+10;
istuck=0;
while (abs(testcw)>1e-4 && i_iter<i_stop  && istuck<4)
    switch DataSource
        case 'BodyWaves'
            [chisqr,J,dv,~,npflg] = BW_calc(Cout);
            sig=sigdat_bw(idfnt);
        case 'SurfaceWaves'
            [chisqr,J,dv,~,npflg] = SW_calc(Cout);
            sig=[sigdat_sw ; dcomp(:)];
        case 'B&SWaves'
           [chisqr,J,dv,~,npflg] = EC_calc(Cout);       
           sig=[sigdat_bw(idfnt);sigdat_sw; dcomp(:)]; 
    end   
    sigm=sig*ones(1,nconst);
    G=J./sigm;  
    alpha=G'*G;  
    a_lam=alpha + diag(diag(alpha)*lam);
    % solve for increments of parameters:
    dC=a_lam\(J'*(dv(:)./sig.^(2)));
    testcw=(oldchisqr-chisqr)/chisqr;
    if (testcw>0 && npflg==0)
        txt=sprintf('%4i       %8.3f      %8.3e       %8.3e     %8.3e',...
            i_iter+InputStrct.opts.funiter,chisqr,testcw,lam,relax);
        if(iter_disp),disp(txt),end
        i_iter=i_iter+1;
        Cbest=Cout;
        Cout=Cout(:)+relax*dC;
        oldchisqr=chisqr;
        lam=lam/10;
        relax=relax*1.25;
        if relax>4,relax=4;end
    elseif (testcw<=0 && npflg==0)
        lam=lam*10; 
        relax=relax/1.25;
        Cout=Cout(:)+relax*dC; 
        istuck=istuck+1;
    elseif (npflg==1)
        lam=lam*10; 
        relax=relax/10; 
        istuck=istuck+1;
    end
    id=find((Cout-lb)<0);
    Cout(id)=lb(id);
    id=find((Cout-ub)>0);
    Cout(id)=ub(id); 
end
Cout=Cbest;
InputStrct.opts.funiter=i_iter;
end

%***********************************************************************************
function [eaout,chisqr]=LM_LSQR_ea(eain,ix,lb,ub)
%function to optimize euler angles using the Levenberg-Marquardt method
% Usage:   eaout=LM_LSQR(eain)
% where eain and eaout are vectors of 3 euler angles for one crystal
% orientation. InputStrct must be an avalable for this routine to work
% calls jacobian_ea
%jmb 3/2010
relax=1e0;
InputStrct.opts.eulerangles(:,ix)=eain;
eabest=eain;
eaout=eain;
testea=1000;
lam=1e-2;
oldchisqr=1e6;
i_iter=InputStrct.opts.funiter;
i_stop=i_iter+10;
istuck=0;
while (abs(testea)>1e-4 && i_iter<i_stop  && istuck<4)
    [chisqr,J,dv,sig]=jacobian_ea(InputStrct,ix,Co);
    sigm=sig*ones(1,3);
    G=J./sigm;  
    alpha=G'*G;  
    a_lam=alpha + diag(diag(alpha)*lam);
    d_ea=a_lam\(J'*(dv(:)./sig.^(2)));
    testea=(oldchisqr-chisqr)/chisqr;
    if (testea>0)
        i_iter=i_iter+1;
        eabest=eaout;
        eaout=eaout(:)+relax*d_ea(:);
        oldchisqr=chisqr;
        lam=lam/5;
        relax=relax*1.5;
        if relax>4,relax=4;end
    elseif (testea<=0)
        lam=lam*5; 
        relax=relax/1.5;
        eaout=eaout(:)+relax*d_ea; 
        istuck=istuck+1;
    end
    id=find((eaout-lb)<0);
    eaout(id)=lb(id);
    id=find((eaout-ub)>0);
    eaout(id)=ub(id);
    InputStrct.opts.eulerangles(:,ix)=eaout;
end
eaout=eabest;
InputStrct.opts.eulerangles(:,ix)=eabest;
InputStrct.opts.funiter=i_iter;
end

%***********************************************************************************
function checkinput
% a set of checks of the input to velocities2Cij
        if (not(isfield(InputStrct,'Data')) && not(isfield(InputStrct,'opts')))
            error('Input structure needs .Data and .opts fields')
        end
        if not(isfield(InputStrct.Data,'dcosflg'))
         InputStrct.Data.dcosflg=0;
        end
        if not(isfield(InputStrct.opts,'pltRange'))
            InputStrct.opts.pltRange=180;
        end
       if not(isfield(InputStrct.opts,'funiter'))
          InputStrct.opts.funiter=0;
       end       
        if (isfield(InputStrct.Data.sample,'SWvelocities'))
            if not(isfield(InputStrct.Data.sample,'BWvelocities'))
                InputStrct.Data.DataSource='SurfaceWaves';
                DataSource='SurfaceWaves';
            else
                InputStrct.Data.DataSource='B&SWaves';
                DataSource='B&SWaves';
            end
        else
            InputStrct.Data.DataSource='BodyWaves';
            DataSource='BodyWaves';
        end
        if not(isfield(InputStrct.Data.sample,'dcos'))
         InputStrct.Data.dcosflg=0;
        end
        sy=InputStrct.Data.sym;
        nc=length(Cin);
        switch sy
            case 'h'
                if nc~=5,error('Hexagonal symmetry requires 5 moduli in input'),end
            case 'o'
                if nc~=9,error('Orthorhombic symmetry requires 9 moduli in input'),end
            case 'm'
                if nc~=13,error('Monoclinic symmetry requires 13 moduli in input'),end
            case 'tri'
                if nc~=21,error('Triclinic symmetry requires 21 moduli in input '),end
            otherwise
                error('symmetry code and number of moduli is not set for this symmetry')
        end
        if not(strcmp(Refine_Cij(1),'y') ||strcmp(Refine_Cij(1),'n') || ...
                strcmp(Refine_Cij(1),'r'))
            error( 'Refine_Cij needs to be y, n, or r')
        end
        [nc,nr]=size(ea);
        ns=InputStrct.Data.nsamp;
        if (InputStrct.Data.dcosflg==0 && (nr~=3 || nc~=ns))
            error('Eulerangle input needs to be 3 x # of samples matrix')
        end
        if not(strcmp(Refine_ea(1),'y') ||strcmp(Refine_ea(1),'n') ||...
                strcmp(Refine_ea(1),'r'))
            error( 'Refine_ea needs to be y, n, or r')
        end
        if not(strcmp(DataSource,'BodyWaves') || strcmp(DataSource,'SurfaceWaves')...
                || strcmp(DataSource,'B&SWaves'))
            error( 'DataSource needs to be ''BodyWaves'' ''SurfaceWaves'' or ''B&SWaves''')
        end
       if not(strcmp(Method,'LM') || strcmp(Method,'NM') || strcmp(Method,'BG'))
            error( 'Method needs to be ''LM'' or ''NM'' or ''BG''')
       end
        
end
end
