function [ Qt , Ut ] = Exp_ST( U0 , time )

  if nargin < 2,  time = 1; end

  d = size( U0 , 1 ); 
  if ~isequal( size(U0) , [d,d] ), error('U0 must be square'); end
  d = d - 1;
  if d < 1, error('invalid dimension'); end
  

  sigma = U0(1,1);        tau = U0(1:d,d+1);
  
  if true
    U0_ = diag( [ sigma*ones(1,d) , 0 ] ); U0_(1:d,d+1) = tau;
    if ~isequal( U0 , U0_ ),  										error('bad structured U0'); end
  end

  
  ntau = sqrt( sum( tau(:).^2  ) );
  

  Qt = zeros( d+1 , d+1 , numel(time) );
  Qk = eye(d+1,d+1); delements = find( Qk ); delements = delements(1:d);

  if nargout > 1
  Ut = zeros( d+1 , d+1 , numel(time) );
  Uk = zeros(d+1,d+1);
  end
  
  if ntau > 0
    
    w   = sqrt( sigma^2 + ntau^2 );
    b   = ( w - sigma  )/ntau;
    denominator = b^2 * exp( 2 * w * time ) + 1;

    for k = 1:numel(time)
      Qk(delements) =  (b^2+1)*   exp(   w * time(k) )       / denominator( k );
      Qk(1:d,d+1)   =   b     * ( exp( 2*w * time(k) ) - 1 ) / denominator( k ) * tau/ntau;

      Qt(:,:,k) = Qk;

      if nargout > 1
      Uk(delements) =      w * (1+b^2) * exp(   w * time(k) ) * ( 2 - denominator( k ) )/denominator( k )^2;
      Uk(1:d,d+1)   =  2 * w * (b+b^3) * exp( 2*w * time(k) ) / denominator( k )^2 * tau/ntau;

      Ut(:,:,k) = Uk / Qk(1);
      end
    end
    
  else
    
    for k = 1:numel(time)
      Qk(delements) = exp( sigma * time(k) );

      Qt(:,:,k) = Qk;

      if nargout > 1
      Uk(delements) = sigma * exp( sigma * time(k) );

      Ut(:,:,k) = Uk / Qk(1);
      end
    end
    
  end

  
end
