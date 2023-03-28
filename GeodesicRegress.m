function [Us,E] = GeodesicRegress( Xs , Y , varargin )
% 

  MAX_ALT_ITS = 40;
  MAX_ITS     = 100;


  N = size(Y,3);

  [varargin,i,EXP_fcn] = parseargs( varargin , 'EXPfcn' , '$DEFS$', [] );
  if ~isa( EXP_fcn , 'function_handle' ), error('''EXPfcn'' must be specified   @(U,varargin)Exp_SxSOxT(U,varargin{:}) '); end

  [varargin,i,LOG_fcn] = parseargs( varargin , 'LOGfcn' , '$DEFS$', [] );
  if ~isa( LOG_fcn , 'function_handle' ), error('''LOGfcn'' must be specified   @(Q)Log(Q)'); end

  [varargin,i,BASIS  ] = parseargs( varargin , 'Basis'  , '$DEFS$', [] );
  if isempty( BASIS ) || ~isnumeric( BASIS ), error('''Basis'' must be specified as a numeric n^2 x k matrix.'); end

  pBASIS = pinv(BASIS);
  k = size( BASIS , 2 );
  n = sqrt( size(BASIS,1) ); In = eye(n,n);

  [varargin,i,W  ] = parseargs( varargin , 'Weight'  , '$DEFS$', 1 );
  W = W(:);
  if isscalar( W )
    W = W*ones(N,1);
  elseif length( W ) ~= N
    error('incorrect size in Weights');
  end


  %preparing the design matrix X
  X = zeros(N,0);
  if ~iscell(Xs), Xs = {Xs}; end
  for c = 1:numel(Xs)
    if isscalar( Xs{c} ),         Xs{c} = Xs{c}*ones(N,1);
    elseif isvector( Xs{c} ),     Xs{c} = Xs{c}(:);
    elseif ndims( Xs{c} ) ~= 2,   error('incorrect input Xs. It must be a matrix!');
    elseif size( Xs{c},1 ) ~= N,  error('incorrect sizes in Xs.');
    end
    X = cat( 2 , X , Xs{c} );
  end
  C = size( X , 2 );


  where = all( isfinite(X) , 2 ) & isfinite( W ) & ~~W;
  X = X(where,:);
  Y = Y(:,:,where);
  W = W(where);
  N = size(Y,3);
  
  iY = nan(size(Y));
  for i = 1:N, iY(:,:,i) = Y(:,:,i) \ In; end

  
  STORED_U = cell(N,C);
  STORED_y = cell(N,C);
  D        = NaN(N,1);

  
  U = zeros(k,C); UU = zeros(k,C);
  for it = 1:MAX_ITS
    Uprev = U;
    
    for c = 1:C
      A = U( : ,     1:(c-1) );
      B = U( : , (c+1):C     );
      U(:,c) = Optimize( @(u)ENER( [ A , u , B ] ) , U(:,c) , 'methods',{'conjugate',ceil(MAX_ALT_ITS/4),'quasinewton',ceil(MAX_ALT_ITS/4)},'ls',{'quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'5'}},'MAX_ITERATIONS',MAX_ALT_ITS),'verbose',0,'noplot');
    end
    
    E = ENER( U );
    fprintf( 'it: %5d   diff in parameters: %16.16g    E: %16.16g\n', it , maxnorm( U , Uprev ) , E );

    if ~mod(it,5) || isequal( U , Uprev )
      U = Optimize( @(U)ENER(U) , U , 'methods',{'conjugate',ceil(MAX_ALT_ITS/4),'quasinewton',ceil(MAX_ALT_ITS/4)},'ls',{'quadratic','golden','quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'5'}},'MAX_ITERATIONS',MAX_ALT_ITS*10),'verbose',0,'plot');
      E = ENER( U );
      fprintf( '*********   diff in parameters: %16.16g    E: %16.16g\n', maxnorm( U , Uprev ) , E );
    end

    if maxnorm( U , Uprev ) < 1e-10, break; end
  end
  
  
  for c = 1:C
    Us{c} = reshape( BASIS*U(:,c) ,[n,n] );
  end
  
  
  function E = ENER( U )
    UU = U;
    for ii = 1:N
      y = In;
      for cc = 1:C
        y = y * getEXP( ii , cc );
      end
      D(ii) = nan;
      try,    D(ii) = LOG_fcn( iY(:,:,ii) * y ); 
      catch,  D(ii) = NaN;
      end
    end
    D = D.^2;
    D( ~isfinite(D) ) = 1e30;
    E = sum( W .* D );
  end

  function yy = getEXP( i , c )
    if isequal( STORED_U{i,c} , UU(:,c) )
      yy = STORED_y{i,c};
    else
      STORED_U{i,c} = UU(:,c);
      yy            = EXP_fcn( reshape( BASIS * ( X(i,c) * UU(:,c) ) ,[n,n] ) );
      STORED_y{i,c} = yy;
    end
  end

  
end
