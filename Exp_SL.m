function [ Qt , Ut ] = Exp_SL( U0 , time )

  if nargin < 2,  time = 1; end

  d = size( U0 , 1 ); 
  if ~isequal( size(U0) , [d,d] ), error('U0 must be square'); end

  if d < 1, error('invalid dimension'); end
  
  if true
    if abs( trace( U0 ) ) > 1e-12,            		error('bad structured V0  ( trace( V0 ) ~= 0 )');         end
  end

  W = U0 - U0.';
  
  Qt = zeros( d , d , numel(time) );
  for k = 1:numel(time)
    Qt(:,:,k) = expm( time(k)*U0.' )*expm( time(k)*W );
  end
    
  if nargout > 1
  Ut = zeros( d , d , numel(time) );
  for k = 1:numel(time)
    Ut(:,:,k) = expm( - time(k)*W ) * U0 * expm( time(k)*W );
  end
  end

end
