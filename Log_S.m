function [d,U] = Log_S( Q )

  d = size( Q , 1 );

  if true
    if ~isequal( size(Q) , [d d] ), error('Q must be square'); end

    
    Q_ = diag(Q);
    if any(Q) <= 0,   error('diag(Q) must be greater than zero.'); end
    
    Q_ = diag(Q_);
    if ~isequal( Q , Q_ ), error('bad structured X0'); end
  end
  
  
  d = norm( log( diag(Q) ) );

  if nargout > 1
    U = diag( log( diag( Q ) ) );
  end
  
end
