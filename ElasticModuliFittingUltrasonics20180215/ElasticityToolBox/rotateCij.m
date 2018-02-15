function cout=rotateCij(cin,atr)  
% USAGE:    cout=rotateCij(cin,atr)
%      where cin == Cijkl or Cij and atr is a rotation matrix
%  if cin is 6x6 matrix, it is converted to 3x3x3x3 for rotation
%
%          J. Michael Brown
%          University of Washington
%          brown@ess.washington.edu             2000
%

%convert to tensor if 6x6 matrix is provided
n=length(size(cin));
if n==2,
       cin=Tnsr2Mtrx(cin);
end

cout=zeros(3,3,3,3);      
for i=1:3,
  for j=1:3,
    for k=1:3,
      for l=1:3,
        for ia=1:3,
          for ja=1:3,
            for ka=1:3,
              for la=1:3,
                 cout(i,j,k,l)=cout(i,j,k,l) +atr(i,ia)*atr(j,ja)*atr(k,ka)*atr(l,la)*cin(ia,ja,ka,la);
              end
            end
          end
        end
      end
    end
  end
end

%convert back to 6x6 matrix  if needed.
if n==2,
    cout=Tnsr2Mtrx(cout);
end



