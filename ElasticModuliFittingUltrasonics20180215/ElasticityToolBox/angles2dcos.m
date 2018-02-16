function   dcos=angles2dcos(a,ea)
% function to determine direction cosines given an orientation matrix
% (crystal coordinates relative to lab coordinates) specified by euler angles "ea" 
% and an angle "a".  The assumption is that measurements are made in lab 
% coordinates in the x-y plane with the angle "a" measured from the x axis.
% Usage: dcos=angles2dcos(a,ea) where a and ea are in degrees.
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2000

 np=length(a);
 atr=eiler(ea);
for j=1:np;
        c=cosd(a(j));s=sind(a(j));
        r(1,1)=c;r(1,2)=s;r(1,3)=0;r(2,1)=-s;r(2,2)=c;r(2,3)=0;r(3,1)=0;r(3,2)=0;r(3,3)=1;
        d=r*atr;
        dcos(j,:)=d(1,:);
end