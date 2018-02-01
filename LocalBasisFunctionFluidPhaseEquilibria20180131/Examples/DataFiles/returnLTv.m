function vout=returnLTv(PT)
   has_parfor = ~isempty(which('parfor'));
    if(iscell(PT))
        P=PT{1};
        T=PT{2};
        flg_grd=1;
        nP=length(P);
        nT=length(T);
        vout=zeros(nP,nT);
    else
        P=PT(:,1);
        T=PT(:,2);
        flg_grd=0;
        nP=length(P);
        nT=length(T);
        vout=zeros(nP,1);
    end
v_guess=[1100 2700];
if(flg_grd)   
    for i=1:nP
        ptmp=P(i);
        if has_parfor
            parfor j=1:nT
              vout(i,j)=fzero(@(v) (ptmp-LTfunc(v,T(j))),v_guess,optimset('disp','off')); 
            end
        else
            for j=1:nT
              vout(i,j)=fzero(@(v) (ptmp-LTfunc(v,T(j))),v_guess,optimset('disp','off')); 
            end
        end
    end
 else
    if has_parfor
        parfor i=1:nP
          vout(i)=fzero(@(v) (P(i)-LTfunc(v,T(i))),v_guess,optimset('disp','off')); 
        end
    else
        for i=1:nP
          vout(i)=fzero(@(v) (P(i)-LTfunc(v,T(i))),v_guess,optimset('disp','off')); 
        end 
    end
        
end

function pc=LTfunc(v,T)
%Lin and Trusler fit coefficients
a=zeros(4,4);
a(1,1)= -6.462689E-01;
a(1,2)= 1.100536E+00;
a(1,3) =1.157852E+00;
a(1,4) =-1.009160E+00;

a(2,1) =-5.483820E-03;
a(2,2)= 1.868520E-02 ;
a(2,3) =-1.823371E-02;
a(2,4) =4.950225E-03;

a(3,1) =2.775202E-06;
a(3,2)= -6.995566E-06;
a(3,3) =3.321712E-06;
a(3,4) =1.298353E-06;

a(4,1) =1.228929E-09;
a(4,2) =-5.279129E-09;
a(4,3) =7.992004E-09;
a(4,4) =-4.125395E-09;
bj=[-4 -1 1 4 6 7];
b= [-2.865737E+02 6.741434E+02 1.475934E+03 -5.363834E+02 2.580254E+02 -8.392639E+01];

X=300*T^-1;
pc=0.0;
Co=sum(b.*X.^-bj);

for j=1:4
   for k=1:4
       pc=pc+a(j,k)*(v-Co).^j*X^(k-1);
   end
end

pc=pc+.1;