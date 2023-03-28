function [ Qt , Ut ] = Exp_GA( U0 , time )

  if nargin < 2,  time = 1; end

  d = size( U0 , 1 ); 
  if ~isequal( size(U0) , [d,d] ), error('U0 must be square'); end
  d = d - 1;
  if d < 1, error('invalid dimension'); end
  
  if true
    if any( abs( U0(d+1,:) ) > 1e-12 ),  										error('bad structured U0'); end
  end

  B = LIEbasis(sprintf('ga(%d)',d));
  m = eye( size(B,2) );
  
  Qt = zeros( d+1 , d+1 , numel(time) );
  Ut = zeros( d+1 , d+1 , numel(time) );
  for k = 1:numel(time)
    [ Qt(:,:,k) , Ut(:,:,k) ] = Exp( Q0 , V0 , B , m , struct('T',time(k)) );
  end

end
