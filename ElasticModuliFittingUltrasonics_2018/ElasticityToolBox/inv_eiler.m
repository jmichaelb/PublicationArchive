function w = inv_eiler(m)
%
% Berechnung der Eulerwinkel aus der Drehmatrix m, wie sie von eiler
% berechnet wird. 
%http://www.physik.uni-osnabrueck.de/resonanz/hjreyher/soljanka/coord_Euler/Euler_in_general/euler_angles_general.htm
%
w=zeros(1,3);
beta=acos(m(3,3));
%----------------------------------- beta > 0 ------------------------------
if abs(sin(beta)) > 1.0e-14
   cos_alpha = m(3,1)/sin(beta);
   cos_alpha = sign(cos_alpha)*min([1 abs(cos_alpha)]);
   sin_alpha = m(3,2)/sin(beta);
   sin_alpha = sign(sin_alpha)*min([1 abs(sin_alpha)]);   
   if sin_alpha >= 0
      if cos_alpha >= 0    
         alpha_grad = acos(cos_alpha)*180/pi;              % 1. Quadrant
      else
         alpha_grad  = acos(cos_alpha)*180/pi;             % 2. Quadrant
      end
   else
      if cos_alpha <= 0.0
         alpha_grad = 360 - acos(cos_alpha)*180/pi;         % 3. Quadrant
      else
         alpha_grad = 360 - acos(cos_alpha)*180/pi;         % 4. Quadrant
      end
   end
   % jetzt dito mit gamma
   cos_gamma  = -m(1,3) / sin(beta);
   cos_gamma = sign(cos_gamma)*min([1 abs(cos_gamma)]);
   sin_gamma  =  m(2,3) / sin(beta);
   sin_gamma = sign(sin_gamma)*min([1 abs(sin_gamma)]);   
   
   if sin_gamma >= 0
      if cos_gamma >= 0    
         gamma_grad = acos(cos_gamma) * 180/pi;
      else
         gamma_grad = acos(cos_gamma) * 180/pi;
      end
   else
      if cos_gamma <= 0.0   
         gamma_grad = 360 - acos(cos_gamma)*180/pi;   
      else
         gamma_grad = 360 - acos(cos_gamma)*180/pi;
      end
   end
else
   % ----- fŸr beta < 1e-10  
   s=1; if  m(3,3) < 0 , s=-1; end
   cos_alphgam = s*m(1,1);    % alphgam steht entweder fŸr alpha + gamma oder alpha - gamma
   cos_alphgam = sign(cos_alphgam)*min([1 abs(cos_alphgam)]);
   sin_alphgam = s*m(1,2);
   sin_alphgam = sign(sin_alphgam)*min([1 abs(sin_alphgam)]);
   if sin_alphgam >= 0
      if cos_alphgam >= 0    
         alpha_grad = acos(cos_alphgam)*180/pi;
      else
         alpha_grad = acos(cos_alphgam)*180/pi;
      end
   else
      if cos_alphgam <= 0
         alpha_grad = 360  - acos(cos_alphgam)*180/pi;   
      else
         alpha_grad = 360  - acos(cos_alphgam)*180/pi;
      end
   end
   gamma_grad = 0;
end

w(1)=alpha_grad; w(2)=beta*180/pi; w(3)=gamma_grad;
