function mask=mk_mask4Gspline(PT)
if(iscell(PT))
    cell_flg=1;
    P=PT{1};
    T=PT{2};
else
    cell_flg=0;
    P=PT(:,1);
    T=PT(:,2);
end
PTmask=[300000 14000
    269000 14000
    237000 14000
    232000 13000
    222000 12400
    213000 10200
    204000 8000
    196000  6500
    180000  4700
    145000 3800
    97000  3600
    85000 3400
    73000 2700
    66000 1900
    62000 1100
    59000  850 % 780
    32000 780%  770
    17000 680 %780 %  760
    6700  490 %560
    5400 450
    2500 320
    1400 280
    1000 245
    300 240
    273 235 %242  %235
    0 242];

if cell_flg
    nP=length(P);
    nT=length(T);
    mask=ones(nP,nT);
    [Pm,Tm]=ndgrid(P,T);
    for i=1:nP
        Tc=interp1(PTmask(:,1),PTmask(:,2),P(i));
        id= T<Tc;
        mask(i,id)=nan;
    end
    id= Pm<8 & Tm>500;
    mask(id)=nan;
    id= Pm<90 & Tm>500 & Tm<900;
    mask(id)=nan;
else
    nT=length(T);
    mask=ones(nT,1);
    for i=1:nT
        Tc=interp1(PTmask(:,1),PTmask(:,2),P(i));
        if(T(i)<Tc),mask(i)=nan;end
        if(P(i)<8 && T(i)>500),mask(i)=nan;end
        if(P(i)<90 && T(i)>500 && T(i)<900),mask(i)=nan;end
    end
end
    
    
    
    