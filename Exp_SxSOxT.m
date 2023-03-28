function [Q,U] = Exp_SxSOxT( U , t )

%   Q = [1,0,0,0,0,0,0;0,1,0,0,0,0,0;0,0,1,0,0,0,0;0,0,0,1,0,0,U(4,7);0,0,0,0,1,0,U(5,7);0,0,0,0,0,1,U(6,7);0,0,0,0,0,0,1];
  if nargin < 2
    t = 1;
  end
  
  Q = eye(7);
  if isa( U , 'dual' ), Q = dual( Q ); end
  Q([46,47,48]) = t*U([46,47,48]);
  
  if isnumeric( U )
    R = expmrot_old( t*U([10,15,2]) );
  else
    R = U(1:3,1:3);
    R([1,5,9]) = 0;
    R = expm( R );
  end
  
  Q(1:3,1:3) = exp( t*U(1) ) * R;

end  
