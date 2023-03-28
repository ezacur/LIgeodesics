function [Y0,V0,E] = GeodesicRegression( Ts , Ys , C , EXP_fun , LOG_fun , BASIS )
if 0
  Ts = (1:5)-3;
  Ys = Exp_ST( [pi -1;0 0] , Ts );
  for n = 1:numel(Ts)
    Ys(:,:,n) = [5 -20;0 1]*Ys(:,:,n);
  end
  
  [Y0,V0,E] = GeodesicRegression( Ys , Ts , [] , @(varargin)Exp_ST(varargin{:}) , @(Q)Log_ST(Q) , LIEbasis('st(1)') );
  U0 = Y0 \ V0;

  
  [Y0,V0,E] = GeodesicRegression( Ys , Ts , { mean(Ts) , KarcherMean( Ys , [] , @(varargin)Exp_ST(varargin{:}) , @(Q)Log_ST(Q) ) } , @(varargin)Exp_ST(varargin{:}) , @(Q)Log_ST(Q) , LIEbasis('st(1)') );
  U0 = Y0 \ V0;
  
  
  for n=1:numel(Ts), disp( Log_ST( ( Y0 * Exp_ST( Ts(n) * U0 ) ) \ Ys(:,:,n) ) ); end
  
end

  if size(Ys,3) ~= numel(Ts), error('incorrect sizes.'); end
  w = ~isfinite(Ts);
  Ts( w )   = [];
  Ys(:,:,w) = [];
  N = numel(Ts);

  pBASIS = pinv(BASIS);
  k = size( BASIS , 2 );
  n = sqrt( size(BASIS,1) ); In = eye(n,n);
  
  

  if  isempty( C )
    T0 = mean( Ts );
    C  = KarcherMean( Ys , [] , EXP_fun , LOG_fun , 'L' );
  elseif iscell( C ) && numel(C) == 2
    T0 = C{1};
    C  = C{2};
  elseif isnumeric(C)
    T0 = C;
  else
    error('no entiendo C   deberia ser { T0 , Y0 }' );
  end
  Ts = Ts - T0;

  yu = nan( 2*k ,1);
  [Y0,Y0] = LOG_fun( C );
  yu(1:k) = pBASIS * Y0(:);

  iC = C \ In;

  U = nan(k,N);
  for i = 1:N
    [Y0,Y0] = LOG_fun( iC * Ys(:,:,i) );
    U(:,i) = pBASIS * Y0(:);
  end

  for c = 1:k
    yu( k+c ) = wlr( Ts , U(c,:).' );
  end
  
  iYs = nan(size(Ys));
  for i = 1:N
    iYs(:,:,i) = Ys(:,:,i) \ In;
  end
  
  D = NaN(N,1);
  Ts = Ts(:).';
  TSuniques = unique( Ts ); TSuniques = TSuniques(:).';
  lastxx = {NaN(k,1),NaN(n,n)};
  
  for it = 1:50
    yuprev = yu;
    yu(     1:k     ) = Optimize( @(y) ENER( y         , yu( (k+1):(2*k) ) ) , yu(     1:k     ) ,'methods',{'conjugate',10,'quasinewton',10},'ls',{'quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'5'}},'MAX_ITERATIONS',20),'verbose',0,'noplot');
    yu( (k+1):(2*k) ) = Optimize( @(u) ENER( yu( 1:k ) , u                 ) , yu( (k+1):(2*k) ) ,'methods',{'conjugate',10,'quasinewton',10},'ls',{'quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'5'}},'MAX_ITERATIONS',20),'verbose',0,'noplot');

    E = ENER( yu(1:k) , yu( (k+1):(2*k) ) );
    fprintf( 'it: %5d   diff in parameters: %0.15g    E: %0.15g\n', it , maxnorm( yu , yuprev ) , E );
    if maxnorm( yu , yuprev ) < 1e-10, break; end
  end
  yu = Optimize( @(yu) ENER( yu(1:k) , yu( (k+1):(2*k) ) ) , yu ,'methods',{'conjugate',50,'quasinewton',50},'ls',{'quadratic','golden','quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'5'}}),'verbose',0,'noplot');
  E = ENER( yu(1:k) , yu( (k+1):(2*k) ) );

  Yc = EXP_fun( reshape( BASIS * yu(     1:k     ) ,[n,n]) );
  Uc =          reshape( BASIS * yu( (k+1):(2*k) ) ,[n,n]);

  if T0 == 0
    Y0 = Yc;
    U0 = Uc;
  elseif T0 > 0
    [Y0,U0] = EXP_fun( -Uc , T0 );
    Y0 = Yc * Y0;
    U0 = -U0;
  elseif T0 < 0
    [Y0,U0] = EXP_fun( Uc , - T0 );
    Y0 = Yc * Y0;
  end
  V0 = Y0 * U0;
  
  
  function E = ENER( h , m )
    if ~isnumeric( h )  ||  ~isequal( h , lastxx{1} )
      lastxx{1} = h;
      lastxx{2} = EXP_fun( reshape( BASIS * h(:) , [n,n]) );
    end
    H = lastxx{2};
    M = reshape( BASIS * m(:) , [n,n] );
    
    for t = TSuniques
      Exp_tM = H * EXP_fun( t * M );
      for i = find( Ts == t )
        try
          D(i) = LOG_fun( iYs(:,:,i) * Exp_tM );
        catch
          D(i) = NaN;
        end
        if ~isfinite( D(i) ), D(i) = 1e30; end
      end
    end
    E = sum( D.^2 );
  end
end
