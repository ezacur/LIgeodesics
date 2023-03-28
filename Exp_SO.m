function [ Qt , Ut ] = Exp_SO( U0 , time )

  if nargin < 2,  time = 1; end

  d = size( U0 , 1 );
  if ~isequal( size(U0) , [d,d] ), error('U0 must be square'); end

  if d < 1, error('invalid dimension'); end
  
  if true
    if maxnorm( U0 + U0.' ) > 1e-10,          		error('bad structured U0  ( U0 ~= skewmatrix )');  end
  end


  Qt = zeros( d , d , numel(time) );
  for k = 1:numel(time)
    Qt(:,:,k) = expm( U0 * time(k) );
  end


  if nargout > 1
  Ut = bsxfun( @times , U0 , ones([1,1,numel(time)]) );
  end

end
