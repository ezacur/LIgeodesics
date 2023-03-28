function [ Qt , Ut ] = Exp_SIM( U0 , time )

  if nargin < 2,  time = 1; end

  d = size( U0 , 1 ); 
  if ~isequal( size(U0) , [d,d] ), error('U0 must be square'); end
  d = d - 1;
  if d < 1, error('invalid dimension'); end
  
  sigma = U0(1,1);                      omega = U0(1:d,1:d); omega(~~eye(d)) = 0;  tau = U0(1:d,d+1);

  if true
    if max(max(abs( omega + omega.' ))) > 1e-12,  error('bad structured U0  ( omega ~= skewmatrix )');        end
    
    U0_ = [ sigma*eye(d) + omega , tau ; zeros(1,d) , 0 ];
    if max(max(abs( U0 - U0_ ))) > 1e-10, 				error('bad structured U0'); end
  end

  dST = [ eye(d)*sigma , tau ; zeros(1,d) , 0 ];
  
  
  Qt = zeros( d+1 , d+1 , numel(time) );
  RR = eye(d+1,d+1);
  if nargout > 1
  Ut = zeros( d+1 , d+1 , numel(time) );
  dRR = zeros(d+1,d+1);
  [SOk,dSOk] = Exp_SO( omega , time );
  [STk,dSTk] = Exp_ST( dST   , time );
  else
  SOk = Exp_SO( omega , time );
  STk = Exp_ST( dST   , time );
  end

  for k = 1:numel(time)
    RR(1:d,1:d) = SOk(:,:,k);
    Qt(:,:,k) = STk(:,:,k) * RR;
    
  if nargout > 1
    dRR(1:d,1:d) = SOk(:,:,k) * dSOk(:,:,k);
    
    Ut(:,:,k) = Qt(:,:,k) \ ( STk(:,:,k) * dSTk(:,:,k) * RR + STk(:,:,k) * dRR );
  end
  end
    
end
