function [d,U] = Log_SxSOxT( Q , rho )

  if nargin < 2, rho = [1;1;1]; end

  t = Q([46,47,48]);
  
  Q = Q(1:3,1:3);
  
  S = realpow( det(Q) , 1/3 );
  s = log(S);

  R = Q/S;
  r = logmrot( R );
  
  d = sqrt( rho(1)*3*s^2 + rho(2)*r(:).'*r(:) + rho(3)*t(:).'*t(:) );
  
  if nargout > 1
    U = zeros(7);
    U([46,47,48]) = t;
    U(1:3,1:3)    = r;
    U([1,9,17])   = s;
  end

end  
