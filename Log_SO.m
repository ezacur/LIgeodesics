function [d,U] = Log_SO( Q )

  d = size( Q , 1 );

  if true
    if ~isequal( size(Q) , [d d] ), error('Q must be square'); end
    
    if maxnorm( Q   * Q.' - eye(d) ) > 1e-8, error('no satisface R*R'' == I'); end
    if maxnorm( Q.' * Q   - eye(d) ) > 1e-8, error('no satisface R''*R == I'); end
  end
  

  U = real( logmrot( Q ) );
  d = sqrt( sum( skewmatrix( U ).^2 ) );

end
