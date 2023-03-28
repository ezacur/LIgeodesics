function [d,U] = Log( Q , B , MetricTensor , Uinitial , varargin )

  n = size(Q,1);
  if isequal( Q , eye(n) )
    d = 0; U = zeros(n); return;
  end

  k = size(B,2);
  if ~isequal( size(Q)            , [ n   , n ] ), error('Q must be square');           end
  if ~isequal( size(B)            , [ n*n , k ] ), error('B must be n*n x k');          end
  if ~isequal( size(MetricTensor) , [ k   , k ] ), error('MetricTensor must be k x k'); end
  
  u = zeros(k,1);  uu = u;
  if nargin == 4, u = B \ Uinitial(:); end

  fEXP = @(u) Exp( reshape( B*u , [n,n] ) , B , MetricTensor , varargin{:} );

  LINSOLVE_opts1 = struct('UT',true,'TRANSA',true);
  LINSOLVE_opts2 = struct('UT',true);
  FACTOR = 0.75;
  
  lastu = u + NaN;
  while ~isequal( lastu  , u )
    [E0,dE,ddE] = ENER(u);
    if E0 < 1e-14,                      break; end

    lastu = u;
    
    G = gaussDirection( dE , ddE );
    
    alpha = 1;
    while 1
      uu = u + alpha * G;
      if isequal( uu , u ), break; end

      E = ENER( uu );
      if E < E0
        
        extraBACKTRAKING = false;
        while 1
          alpha = alpha * FACTOR;
          un = u + alpha * G;
          En = ENER( un );
          if En < E
            uu = un;
            E  = En;
          else, break;
          end
          extraBACKTRAKING = true;
        end

        if ~extraBACKTRAKING, break; end
        
        alpha = alpha / FACTOR;
        while 1
          alpha = alpha / FACTOR;
          un = u + alpha * G;
          En = ENER( un );
          if En < E
            uu = un;
            E  = En;
          else, break;
          end
        end
        
        break;
      else
        alpha = alpha * 0.5;
      end
    end

%     fprintf('E: %20.20g\n',E);
    u  = uu;
    if E < 1e-14, break; end
  end

  U = reshape( B*u , [n,n] );
  d = sqrt( u.' * MetricTensor * u );

  
  function [E,dE,ddE] = ENER( u )
    X = fEXP( u );
    R = X(:) - Q(:);

    if maxnorm( R./Q(:) ) < 1e-8
      E = -1; dE = []; ddE = [];
      return;
    end
    E = R.'*R;
    
    if nargout > 1
      D  = NumericalDiff( @(u)fEXP(u) , u , 'c' );
      dE = 2*R.'*D;
    end
    if nargout > 2
      ddE = 2*D.'*D;
    end
  end  

  function G = gaussDirection( dE , ddE )
    G = - dE(:);
    try
      L = chol(ddE);
      Ld = diag(L);
      L( ~~eye(size(L)) ) = max( Ld , max( Ld ) * 2.3e-16 );
      
      G = - linsolve( L , linsolve( L , dE(:) , LINSOLVE_opts1 ) , LINSOLVE_opts2 );
    catch
      [evecs,evals] = eig( ddE );
      evals = diag(evals);
      if max( evals ) <= 0, return; end
      evals = max( evals , max(evals)*1e-8 );
      G = - ( evecs * ( ( evecs.' * dE(:) ) ./ evals ) );
    end
  end

end
