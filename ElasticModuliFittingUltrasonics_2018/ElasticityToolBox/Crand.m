function c=Crand(TrustRegion)
% function to create random elastic constants within the trust region
%Usage: c=Crand(TrustRegion)
%       where TrustRegion is a # of constants by 2 matrix that defines the
%       range of acceptable values of the constants
% The function checks that the resulting constants are positive definite
% calls makeCij,  only orthorhombic, monoclinic and triclinic symmetry at
% this point.  Easy to fix
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013


lb=TrustRegion(:,1);
db=diff(TrustRegion(:,1:2)')';

nc=length(lb);
testnp=-1;

while testnp<0
c=lb+rand(nc,1).*db;
if nc==9,
    cij=Ci2Cij(c,'o');
elseif nc==13
    cij=Ci2Cij(c,'m');
elseif nc==21
    cij=Ci2Cij(c,'tri')';
end
testnp=min(eig(cij));
end
    
