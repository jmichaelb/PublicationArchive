function c=Tnsr2Mtrx(c_in)
% function to convert elastic moduli between tensor C(i,j,k,l) 
%    and matrix C(i,j)
% Usage:  c_out=Tnsr2Mtrx(c_in)
%            if c_in is 6x6 matrix c_out is 3x3x3x3 array
%            if c-in is 3x3x3x3 array  c-out is 6x6
%  
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             7/2013

length(size(c_in));
if(length(size(c_in))==2),
  c=zeros(3,3,3,3);
  for i=1:3,
    for j=1:3,
        for k=1:3,
            for l=1:3,
                m=9-i-j;
                if (i==j), m=i; end
                n=9-k-l;
                if (k==l), n=k; end
                c(i,j,k,l)=c_in(m,n);
            end                  
        end
     end
   end
else 
    c=zeros(6,6);
    for i=1:3,
        for j=1:3,
            for k=1:3,
                for l=1:3,
                     m=9-i-j;
                     if (i==j), m=i; end
                     n=9-k-l;
                     if (k==l), n=k; end
                     c(m,n)=c_in(i,j,k,l);
                end
            end
        end
    end
end


end  % of T2M

    