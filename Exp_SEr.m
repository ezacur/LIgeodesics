function [ Qt , Ut ] = Exp_SEr( U0 , time )

  if nargin < 2,  time = 1; end
  
  d = size( U0 , 1 ); 
  if ~isequal( size(U0) , [d,d] ), error('U0 must be square'); end
  d = d - 1;
  if d < 1, error('invalid dimension'); end
  
  omega = U0(1:d,1:d);         tau = U0(1:d,d+1);
 
  if true
    if max(max(abs( omega + omega.' ))) > 1e-12,  error('bad structured U0  ( U0 ~= skewmatrix )'); end

    U0_ = [ omega , tau ; zeros(1,d) , 0 ];
    if ~isequal( U0 , U0_ ), 											error('bad structured U0'); end
  end

  
  Qt = zeros( d+1 , d+1 , numel(time) );
  Qk = eye( d+1 , d+1 );
  for k = 1:numel(time)
    Qk(1:d,1:d) = expm( omega * time(k) );
    Qk(1:d,d+1) = Qk(1:d,1:d) * tau * time(k);
    
    Qt(:,:,k)   = Qk;
  end


  if nargout > 1
    error('not implemented');
  Ut = bsxfun( @times , U0 , ones([1,1,numel(time)]) );
  for k = 1:numel(time)
    Ut(1:d,d+1,k) =  Qt(1:d,1:d,k) \ tau;
%     Ut(1:d,d+1,k) =  expm( - omega * time(k) ) * tau;
  end
  end

end
