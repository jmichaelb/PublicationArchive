function Dreh_M = eiler(aa);
% Dreh_M=eiler(a,b,c), liefert Drehmatrix für Drehung um Eulersche Winkel
% 1) um a um z, 2) um b um y' und 3) um c um z'' (Whittaker, Tinkham und Grachew)
% vektor(im gedrehten System)=Dreh_m*vektor(im Originalsystem)
% gedrehter vektor(im alten System)=Dreh_m' * Originalvektor
% In den Zeilen von Dreh_m stehen die Komponenenten (im Ursprungssystem)
% der Achsen des gedrehten Systems
% http://www.physik.uni-osnabrueck.de/resonanz/hjreyher/soljanka/coord_Euler/Euler_in_general/euler_angles_general.htm

Dreh_M=zeros(3);
a=aa(1)*pi/180;
b=aa(2)*pi/180;
c=aa(3)*pi/180;
Dreh_M(1,1)=cos(a)*cos(b)*cos(c)-sin(a)*sin(c);
Dreh_M(1,2)=sin(a)*cos(b)*cos(c)+cos(a)*sin(c);
Dreh_M(1,3)=-sin(b)*cos(c);
Dreh_M(2,1)=-cos(a)*cos(b)*sin(c)-sin(a)*cos(c);
Dreh_M(2,2)=-sin(a)*cos(b)*sin(c)+cos(a)*cos(c);
Dreh_M(2,3)=sin(b)*sin(c);
Dreh_M(3,1)=cos(a)*sin(b);
Dreh_M(3,2)=sin(a)*sin(b);
Dreh_M(3,3)=cos(b);