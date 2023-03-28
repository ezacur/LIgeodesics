function [ Q1 , U1 ] = Expr( U0 , B , m , varargin )

  Exp_left_FCN = []; T = [];
  if      nargin == 2 &&  isa( B ,'function_handle' )
    Exp_left_FCN = B;
  elseif  nargin == 3 &&  isa( m ,'function_handle' )
    Exp_left_FCN = m;
    T  = B;
  elseif nargin < 3
    error('use Expr( U0 , Basis , MetricTensor , varargin )');
  end
  
  if ~isempty( Exp_left_FCN )

    if isempty( T ), T = 1; end
    if ~isvector( T ), error('T must be a vector'); end
    
    [Q1,U1] = feval( Exp_left_FCN , -U0 , T );

  else
    
    [Q1,U1] = Exp( -U0 , B , m , varargin{:} );
    
  end

  U1 = - U1;
  for k = 1:size(Q1,3)
    Q1(:,:,k) = Q1(:,:,k) \ eye( size(Q1,1)*[1 1] );
  end

end
