function [d,U] = Log_SEr( Q )

  d = size( Q , 1 );
  d = d-1;

  if true
    if ~isequal( size(Q) , [d+1 d+1] ), error('Q must be square'); end
    
    if any( Q(d+1,1:d) ),   error('la ultima fila no es zero'); end
    if Q(d+1,d+1) ~= 1,     error('el ultimo elemento no es 1'); end
    
    if maxnorm( Q(1:d,1:d)   * Q(1:d,1:d).' - eye(d) ) > 1e-10, error('no satisface R*R'' == I'); end
    if maxnorm( Q(1:d,1:d).' * Q(1:d,1:d)   - eye(d) ) > 1e-10, error('no satisface R''*R == I'); end
  end
  
  U = zeros( d+1 , d+1 );
  U(1:d,1:d) = logmrot( Q(1:d,1:d) );
  U(1:d,d+1) = Q(1:d,1:d) \ Q(1:d,d+1);
  
  d = sqrt( sum( skewmatrix( U(1:d,1:d) ).^2 ) + sum( U(1:d,d+1).^2 ) );
    
end
