function [ Qt , Ut ] = Exp_S( U0 , time )

  if nargin < 2,  time = 1; end

  d = size( U0 , 1 ); 
  if ~isequal( size(U0) , [d,d] ), error('U0 must be square'); end

  if d < 1, error('invalid dimension'); end
  
  sigma = diag( U0 );
  
  if true
    if ~isequal( U0 , diag( sigma ) ), 											error('bad structured U0'); end
  end

  Qt = zeros( d , d , numel(time) );
  for k = 1:numel(time)
    Qt(:,:,k) = diag( exp( sigma * time(k) ) );
  end

  if nargout > 1
  Ut = bsxfun( @plus , U0 , zeros([1,1,numel(time)]) );
  end

end
