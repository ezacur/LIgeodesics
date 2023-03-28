function [ Qt , Ut ] = Exp_T( U0 , time )

  if nargin < 2,  time = 1; end

  d = size( U0 , 1 ); 
  if ~isequal( size(U0) , [d,d] ), error('Q0 must be square'); end
  d = d - 1;
  if d < 1, error('invalid dimension'); end
  
  tau = U0(1:d,d+1);
  
  if true
    U0_ = zeros( d+1,d+1 ); U0_(1:d,d+1) = tau;
    if ~isequal( U0 , U0_ ), 											error('bad structured U0'); end
  end

  Qt = bsxfun( @times , U0 , vec( time , 3 ) ) + bsxfun( @plus , eye(d+1) , zeros([1,1,numel(time)]) );

  if nargout > 1
  Ut = bsxfun( @plus , U0 , zeros([1,1,numel(time)]) );
  end

end
