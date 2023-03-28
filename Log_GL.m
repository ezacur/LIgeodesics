function [d,U] = Log_GL( Q , INITIAL )

if 0
%%  

for init = 0:10
  t=nan(100,11); e=t;
  randn('state',0);
  for kk=1:numel(t)
    U = randn(2);
    U = 0.1*U/norm(U,'fro');
    Q = Exp_GL(U);
    tic;
    for i=1:4
      [u,U]=Log_GL(Q,{init,nan});
    end;
    tt = toc;
    t(kk,init+1) = tt/i;
    e(kk,init+1) = maxnorm(Q-Exp_GL(U));
  end

  fprintf('init:%3d  , mean(t) = %.15g   median(t) = %.15g    range(t) = [ %.15g , %.15g ]\n',init, mean(t),median(t),min(t),max(t) )
end

%%
end

if 0
%%
  results = NaN(500,2);
  for it = 1:size(results,1)
    randn('state',it);
    disp(it)
    U0 = randn(3)*1;
    Q  = Exp_GL( U0 );

    tic;
    for tt=1:15
      U = Log_GL( Q ); %disp( maxnorm( U , U0 ) )
    end
    t = toc; t = t/tt;

    results(it,:) = [ t , log10( maxnorm( U , U0 ) ) ];
  end

end

  n  = size(Q,1);
  if ~isequal( size(Q) , [n n] ), error('Q must be square'); end
  if det(Q) <= 0                , error('det(Q) must be positive'); end
  
  n2 = n*n;
  n3 = n2*n;

  Kn      = comm(n);
%   Knt     = Kn.';
  In      = eye(n);
  Inn     = eye(n*n);
  Inn_Kn  = Inn-Kn;
%   Inn_Knt = Inn_Kn.';
  
%   kron1   = kron( In , Kn , In )*kron( Inn , In(:) );
%   kron2   = kron( In , Kn , In )*kron( In(:) , Inn );

  LINSOLVE_opts1 = struct('UT',true,'TRANSA',true);
  LINSOLVE_opts2 = struct('UT',true);
%   LINSOLVE_opts3 = struct('SYM',true,'POSDEF',true);
  FACTOR = 0.75;
  
  
  if nargin < 2, INITIAL = {}; end
  if ~iscell(INITIAL)
    if isvector( INITIAL )
      INITIAL = num2cell( INITIAL );
    else
      INITIAL = {INITIAL};
    end
  end
  INITIAL = [ INITIAL(:).' , 1 , 0 , num2cell( 2:20 ) ];
  for i = 2:numel(INITIAL)
    if any( arrayfun(@(j)isequal(INITIAL{j},INITIAL{i}), 1:i-1 ) )
      INITIAL{i} = [];
    end
  end
  INITIAL = INITIAL( ~cellfun('isempty',INITIAL) );


  
  [R0_U,R0_D,R0_V] = svd( Q );
  R0 = R0_U * R0_V.';% * realpow( detQ , 1/n );
  function U = getINITIAL()
    if isempty( INITIAL )
      U = randn(n);
%       fprintf( 'random INITIAL: %s\n',uneval(U) );
    elseif ~isscalar( INITIAL{1} )
      U = INITIAL{1}; INITIAL(1) = [];
%       fprintf( 'INITIAL: %s\n',uneval(U) );
    elseif isscalar( INITIAL{1} ) && isnan( INITIAL{1} )
%       fprintf( 'no more initials' );
      error('no more initials');
    else
%       fprintf( 'INITIAL case %d: ', INITIAL{1} );
      switch INITIAL{1}
        case 0
          U = zeros(n);
        case 1
          U = real( logm( R0 * trace( R0_D )/n ) );
        case 2
          U = real( logm( R0 ) + In*log( det(Q) )/n );
        case 3
          U = real( logm( R0 )/2 + In*log( det(Q) )/n );
        case 4
          U = real( logm( diag(diag( Q * R0.' )) ) );
          U = real( (U+U.') + logm(R0) )/2;
        case 5
          U = real( logm( Q * R0.' ) );
          U = real( (U+U.') + logm(R0) )/2;
        case 6
          U = logm( R0_U*R0_D*R0_U.' ) + logm( R0_U*R0_V.' )/2;
        case 7
          U = real( logm( Q ) );
          U = real( logm( Q*expm( U.' - U ) ) ).';
        case 8
          U = real( logm(Q) );
        case 9
          U = 0.9*Q;
        case 10
          U = In;
        case 11
          U = Q - In;
        case 12
          [R0,T0] = qr(Q.');
          if det(T0) < 0, T0(1,:) = -T0(1,:); end
          [a,b,c] = svd( T0 );
          U = logm( c * diag([ ones(size(T0,1)-1,1) ; det(c*a.')]) * a.' * T0 );

        otherwise
%           fprintf( 'random matrix: ' );

          U = randn(n);
      end
      INITIAL(1) = [];
%       fprintf( 'U: %s\n', uneval(U) );

    end
    U = U(:);
  end
  

  try,   U = getINITIAL();
  catch, error('INITIAL not specified'); end
  

  lastU = U + NaN;
  while ~isequal( lastU , U )
    [E0,dE,ddE] = ENER(U);
    if E0 < 1e-14
      CL = checkConjugatePoints( U );
      if CL < 1
        U = CL/2*U;
        continue;
      else
        break;
      end
    end
    
    lastU = U;
    
    G = gaussDirection( dE , ddE );

    alpha = 1;
    while 1
      UU = U + alpha * G;
      if isequal( UU , U )
        try
          UU = getINITIAL();
        catch
          d = sqrt( U.' * U );
          warning('not congerved');
          return;
        end

        E0 = ENER(UU);
        break;
      end

      E = ENER( UU );
      if E < E0
        
        extraBACKTRAKING = false;
        while 1
          alpha = alpha * FACTOR;
          Un = U + alpha * G;
          En = ENER( Un );
          if En < E
            UU = Un;
            E  = En;
          else, break;
          end
          extraBACKTRAKING = true;
        end

        if ~extraBACKTRAKING, break; end

        alpha = alpha / FACTOR;
        while 1
          alpha = alpha / FACTOR;
          Un = U + alpha * G;
          En = ENER( Un );
          if En < E
            UU = Un;
            E  = En;
          else, break;
          end
        end
        
        break;
      else
        alpha = alpha * 0.5;
      end
    end

    U = UU;

%     fprintf('E: %20.20g\n',E);
    if E < 1e-10
%       CL = checkConjugatePoints( U );
%       if CL < 1
%        U = CL/2*U;
%        continue;
%       else
        break;
%       end
    end

  end

  U = reshape( U , [n,n] );
  d = sqrt( U(:).'*U(:) );

  
  function [E,dE,ddE] = ENER( U )
    U = reshape( U , [n,n]);
    if nargout == 1

      X = ExpGL( U );
      R = X(:) - Q(:);

      if max(abs( nonans( R ./ Q(:) ) )) < 1e-10
        E = -1; dE = []; ddE = []; return;
      end
      
    else
      
      [X,J] = ExpGL( U );
      R = X(:) - Q(:);

      if max(abs( nonans( R ./ Q(:) ) )) < 1e-10
        E = -1; dE = []; ddE = []; return;
      end

      dE  = 2*R.'*J;
      ddE = 2*J.'*J;
      
    end

    E = R.'*R;
  end
  
  function [X,J] = ExpGL( U )
    W       = U - U.';
    eUt     = expm( U  ).';
    eW      = expm( W  );

    X = eUt*eW;
    
    if nargout > 1
      deUt    = d_expm( U.' );
      deW     = d_expm( W   );

      J = reshape( reshape( deUt*Kn , [n3,n] )*eW ,[n2,n2] ).' + reshape( eUt*reshape( deW*Inn_Kn , [n,n3]) , [n2,n2] );
    end
  end

  function [ detJ , e ] = determinantJacobian( U )
    U = reshape( U , [n,n] );
    
    [X,J] = ExpGL(U);
    detJ = det( J );
    
    if nargout > 1
      R = X(:) - Q(:);
      e = R(:).'*R(:);
    end
  end

  function t = checkConjugatePoints( U )
    [t,v] = Optimize1D( @(t)determinantJacobian( (t-1)*U )^2 , [1 2],[],'methods','exhaustive','min',1,'max',2 );
    t = t-1;
    if v > 0, t = 2; end
        
    if 0
      %%    
      ts = linspace(-0.01,2.1,250); detsJ = ts*0; eners = detsJ;
      for t = 1:numel(ts)
        [detsJ(t),eners(t)] = determinantJacobian( ts(t)*U );
      end
      figure;
      subplot(211);
      plot( ts , detsJ , '.-r' ); hline(0);
      subplot(212);
      plot( ts , eners , '.-b' ); hline(0);
      %%
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