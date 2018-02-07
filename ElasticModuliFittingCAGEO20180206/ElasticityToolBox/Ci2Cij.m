function Cout=Ci2Cij(Cin,sym)
% function to convert between vector and matrix representation of elastic
% moduli.  Usage:
%        Cout=Ci2Cij(Cin,sym)
% where Cin (as vector) are in the (cyclic) order 
%   isotropic            K and G
%   Cubic................C11 C12 C44 
%   Hexagonal............C11 C12 C13 C33 C44 
%   Trigonal(6 Cij's)....C11 C12 C13 C14 C33 C44 
%   Tetragonal(6 Cij's)..C11 C12 C13 C33 C44 C66 
%   Orthorhombic.........C11 C12 C13 C22 C23 C33 C44 C55 C66
%   Monoclinic(MO2)......C11 C12 C13 C15 C22 C23 C25 C33 C35 C44 C46 C55 C66
%   Triclinic(21 Cij's)..C11,C12,C13,...,C22,C23,...C33,C34,...,C44 etc.%
%   symmtry codes are 'i' 'c' 'h' 'tr6' 'tet6' 'o' 'm' 'tri'
%   returns vector if input is 6x6 matrix
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013    
%   
[n,m]=size(Cin);
C=zeros(6,6);
nc=max([n m]);
flg=0;
if min([n m])==1
    flg=1;
else
    if (n~=6 && m~=6)
        error('Matrix of moduli must be 6x6')
    end
end

switch sym
    case 'i' %  Isotropic
       pos=[1  2  3  8  9 15 22  29 36
            1  2  2  1  2  1  3   3  3];
        if flg
            if nc~=2,error('wrong number of moduli for isotropic symmetry'),end
            Ct=[Cin(1)+4/3*Cin(2) Cin(1)-2/3*Cin(2)*1.00001 Cin(2)];
            C(pos(1,:))=Ct(pos(2,:));
            Cout=C+tril(C,-1)';
        else
            Cout=[Cin(1)-4/3*Cin(36) Cin(36)];
        end
    case 'c'  % cubic  
        pos=[ 1  2  3 8  9 15 22 29 36
              1  2  2 1  2  1  3  3  3];
        if flg
         if nc~=3,error('wrong number of moduli for cubic symmetry'),end
         C(pos(1,:))=Cin(pos(2,:));
         Cout=C+tril(C,-1)';
        else
          Cout=zeros(3,1);
          Cout=Cin([1 2 22]);
        end
    case 'h' % hexagonal 5 constants */
        pos=[ 1  2  3  8  9  15  22  29  
              1  2  3  1  3   4   5   5];
        if flg
            if nc~=5,error('wrong number of moduli for hexagonal symmetry'),end
            C(pos(1,:))=Cin(pos(2,:));
            C(36)=(Cin(1)-Cin(2))/2;
            Cout=C+tril(C,-1)' ;
        else
            Cout=zeros(1,5);
            Cout(pos(2,:))=Cin(pos(1,:));
        end
    case 'tr6' % trigonal 6 constants
        pos=[ 1 2 3  4  8 9 15 22 29  
              1 2 3  4  1 3  5  6  6 ];
        if flg
            if nc~=6,error('wrong number of moduli for trigonal6 symmetry'),end
            C(pos(1,:))=Cin(pos(2,:));
            C(10)=-C(4);
            C(36)=(Cin(1)-Cin(2))/2;
            C(30)=C(4)/2;
            Cout=C+tril(C,-1)' ;
        else
            Cout=zeros(1,6);
            Cout(pos(2,:))=Cin(pos(1,:));
        end
    case 'tet6'  % tetragonal 6 constants    
        pos=[1 2 3 8 9 15 22 29 36
             1 2 3 1 3  4  5  5  6];
         if flg
            if nc~=6,error('wrong number of moduli for tetragonal6 symmetry'),end
            C(pos(1,:))=Cin(pos(2,:));
            Cout=C+tril(C,-1)' ;
         else
             Cout=zeros(1,6);
             Cout(pos(2,:))=Cin(pos(1,:));
         end
    case 'o'  % orthorhombic 9 constants
         pos=[1 2 3 8 9 15 22 29 36]; %
         if flg
            if nc~=9,error('wrong number of moduli for orthorhombic symmetry'),end
            C(pos)=Cin;
            Cout=C+tril(C,-1)';
         else
             Cout=zeros(1,9);
             Cout=Cin(pos);
         end
    case 'm' 
         pos=[1 2 3 5 8 9 11 15 17 22 24 29 36
              1 2 3 4 5 6  7  8  9 10 11 12 13]; 
          if flg
              if nc~=13,error('wrong number of moduli for monoclinic symmetry'),end
              C(pos(1,:))=Cin(pos(2,:));
              Cout=C+tril(C,-1)';    
          else
              Cout=zeros(1,13);
              Cout(pos(2,:))=Cin(pos(1,:));
          end
    case 'tri'  % triclinic 21 constants        
         pos=[1:6 8:12 15:18 22:24 29:30 36];
          if flg
              if nc~=21,error('wrong number of moduli for triclinic symmetry'),end
              C(pos(1,:))=Cin;
              Cout=C+tril(C,-1)';
          else
              Cout=zeros(1,21);
              Cout=Cin(pos(1,:));
          end
    otherwise
        error('Unsupported Symmetry for function makeCij')
end

              

        
  