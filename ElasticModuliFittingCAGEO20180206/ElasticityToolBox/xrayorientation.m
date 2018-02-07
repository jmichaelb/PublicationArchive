function [OM,OA,ce,angles,volume]=xrayorientation(U,T,EA,offset,aflg)
% Usage:
%  [OM,OA,edges,angles,volume]=xrayorientation(U,T,euler,offset,atr)
% returns orientation and cell parameter data - in ISS lab coordinates
% OM has cartesian coordinates centered on crystal axes a*,b, and c*
% OA has actual crystal axes
% ce are lengths of crystal axes
% angles are angles between axes
% volume is the cell volume
% U and T are the matrices in the RMAT file.  U has reciprocal vectors and T is the transformation
% in real space between the primitive and conventional cell
% euler = [phi omega chi] are the angles necessary to rotate sample into a known orientation relative to
% the laboratory coordinate system
% offset is the rotation angle needed to register the x-ray lab coordinates to the ISS coordinates.
% aflg is flag for flipping slide for double-sided polish.

if nargin==4,
    atr=eye(3,3);
else 
    if aflg=='y',
        atr=[-1 0 0;0 1 0;0 0 -1]; % flip slide over by rotating about y (vertical) axis
    else
        atr=eye(3,3);
    end
end

%  transform to real space primitive lattice
     Mp=inv(U');
%  convert to conventional lattice
     Mc=Mp*T';  
    volume=abs(det(Mc));
     
%  calculate reciprocal conventional cell and metric tensor
     Uc=inv(Mc');
     G= inv(Uc'*Uc);
     ce=sqrt(diag(G));

%angles between conventional axes:
    beta=acos(G(1,3)/ce(1)/ce(3))*180/pi;
    alpha=acos(G(2,3)/ce(2)/ce(3))*180/pi;
    gamma=acos(G(1,2)/ce(1)/ce(2))*180/pi;   
    angles=[alpha beta gamma];

% direct lattice vectors in lab frame, placed as columns in matrix
% rotate a specified crystal direction into a known lab direction
    OA=rotateNonius(EA(1),EA(2),EA(3),Mc) ;

%  x-ray work places slide in z-x plane - our calculations use slide in the
%  x-y plane -> mapping x->x  z->y  y->-z
    r=[1  0  0
       0  0  1
       0 -1  0];
    OA=r*OA;
% a matrix of unit vectors parallel to the crystal axes
    O=[OA(:,1)/norm(OA(:,1)) OA(:,2)/norm(OA(:,2)) OA(:,3)/norm(OA(:,3))];
% now to convert crystal axes to cartesian coordinates "CC"
% set bCC parallell to b
    bCC=O(:,2);
% find a* 
    aCC=cross(O(:,2),O(:,3));
% find c
    cCC=cross(aCC,bCC);
% and the orientation matrix is:
    OM=[aCC(:)/norm(aCC) bCC(:)/norm(bCC) cCC(:)/norm(cCC)];  
    OM=atr*OM;
% now rotate orientation matrix into the ISS lab frame
%  positive angle of rotation moves slide XY clockwise, looking down Z
    R=offset*pi/180;
    XYZ_in_labframe=[ [cos(R) sin(R) 0]' [-sin(R) cos(R) 0]' [0 0 1]' ]; 
    OM=XYZ_in_labframe*OM;
    OA=XYZ_in_labframe*OA;
    
% Nested function    **********************************************************************************
 
function rotated_matrix=rotateNonius(phi,chi,omega,zero_matrix)
% function rotated_matrix=rotateNonius(phi,chi,omega,zero_matrix)
% takes the orientation matrix given by nonius program "xmatrix", stored in a .rmat file
% e.g. xmatrix MgO.rmat nostandardize notransformtrigonal 
% and finds the corresponding matrix after rotation of the crystal by "Euler" angles phi, chi, omega
% where phi is around the laboratory Z axis in the -x*y direction, chi is about the laboratory X axis in the
% -y*z direction, and omega is a final rotation about the lab Z. 
% the .rmat, orientation matrix has columns corresponding to inverse lattice vectors a*, b* and c*
% and rows corresponding to laboratory X, Y and Z.

phi=phi*pi/180;
chi=chi*pi/180;
omega=omega*pi/180;

phirotate=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1];

chirotate=[1 0 0; 0 cos(chi) sin(chi); 0 -sin(chi) cos(chi)];

omegarotate=[cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

rotated_matrix=omegarotate*chirotate*phirotate*zero_matrix;

% from Bob Downs notes for 4 circle - not right for this situation
% co=cos(omega);so=sin(omega),cc=cos(chi);sc=sin(chi);cp=cos(phi);sp=sin(phi);
% 
% rotated_matrix=[co*cc*cp-so*sp co*cc*sp+so*cp -co*sc;-so*cc*cp-co*sp -so*cc*sp+co*cp so*sc;sc*co sc*sp co]*zero_matrix;
 
