function  out=KG_calc(c,Mc,Ms,sym)
% this propagates errors using the standard linear systems.
%  given vector of data x and its covariance Mx and given a quantity that
%  is a linear combination of the x's F (F=A'*x), then the covariance in F
%  is given as Mf=A'*Mx*A and the uncertainties are sqrt(diag(Mx)).
% Here the problem is uncertainties for K and G calculated from Cij and its
% covariance.  This uses the Voigt Reuss average and assumes a vector of Cij in
% cyclic order.
% Usage:   out=KG_calc(c,Mc,Ms,sym)  where c is  a vector of
% constants and Mc and Ms are the covariance matrixes.  out has K and G in V and Ruess limits with
% uncertainties 
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

f1=1/9;
f2=1/15;

Av= [1  1  % 11 1
    1  1  % 22  2
    1  1  % 33  3
    0  3  % 44  4
    0  3  % 55  5
    0  3  % 66  6
    2 -1  % 12  7
    2 -1  % 13  8
    0  0  % 15  9
    2 -1  % 23 10
    0  0];
Av(:,1)=Av(:,1)*f1;
Av(:,2)=Av(:,2)*f2;

Ar=[1  4
    1  4
    1  4
    0  3
    0  3
    0  3
    2 -4
    2 -4
    0  0
    2 -4
    0  0];
Ar(:,2)=Ar(:,2)/15;

switch sym
    case 'c'
        pos=[1 7 4];
    case 'h'
        c=[c (c(1)-c(2))/2];
        pos=[1 7 8 1 10 3 5 5 6];
    case 'o'
        pos=[1 7 8 2 10 3 4 5 6 ];
    case 'm'
        pos=[1 7 8 9 2 10 11 3 11 4 11 5 6];
    case 'tri'
        pos=[1 7 8 11 11 11 2 10 11 11 11 3 11 11 11 4 11 11 5 11 6];
    otherwise
        display('unsupported symmetry in KG_calc')
end

temp=Av(pos,:)'*c(:);
Kv=temp(1);Gv=temp(2);
temp=2*sqrt(diag(Av(pos,:)'*Mc*Av(pos,:)));  % always work in 2 sigma values
dKv=temp(1); dGv=temp(2);

% need to make same calculation for Reuss average
 s=Ci2Cij(inv(Ci2Cij(c,sym)),sym);
 temp=Ar(pos,:)'*s(:);
 Kr=1/temp(1);Gr=1/temp(2);
 temp=2*sqrt(diag(Ar(pos,:)'*Ms*Ar(pos,:)));
 dSr=temp(1);dGr=temp(2); 
 out= [Kv dKv Kr dSr*Kr^2 Gv dGv Gr dGr*Gr^2] ;
 out=round(out*10)/10;
