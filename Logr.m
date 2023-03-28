function [d,U] = Logr( Q , varargin )

  Q = Q \ eye(size(Q));
  [d,U] = Log( Q , varargin{:} );
  U = -U;

end
