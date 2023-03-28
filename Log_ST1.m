function [dist,U] = Log_ST1( Q , M )

  if nargin < 1, error('at least 1 argument is expected'); end
  if nargin < 2 || isempty( M ), M    = eye(2); end

  if max(abs( vec( M - M.' ) ) ) > 1e-10
    error('M must be symmetric');
  end

  

  s =  sqrt( det(M) )/M(2,2);
  t =  M(1,2)/M(2,2);
  
  Y = [s , t;0 , 1]; iY = [1/s , -t/s ; 0 , 1];
  
  Q = iY * Q * Y;
  
  [dist,U] = Log_ST( Q );

  U = Y * U * iY;
  
  
  dist = sqrt( U(:).' * [1,0,0,0;0,0,1,0].' * M * [1,0,0,0;0,0,1,0] * U(:) );
  
end
