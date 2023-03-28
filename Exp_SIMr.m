function [ Qt , Ut ] = Exp_SIMr( U0 , varargin )

  [Qt,Ut] = Exp_SIM( -U0 , varargin{:} );
  
  I = eye(size(Qt(:,:,1)));
  for t = 1:size(Qt,3)
    Qt(:,:,t) =   Qt(:,:,t)\I;
    Ut(:,:,t) = - Ut(:,:,t);
  end

end
