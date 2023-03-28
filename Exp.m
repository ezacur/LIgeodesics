function [Qt,Ut,ITERATIONS] = Exp( U0 , B , MetricTensor , OPTS )
%Exp   :   Riemannian exponential under a left-invariant metric.
%
%   Q1 = Exp( U0 , matBases , MetricTensor )
%     Compute the Riemannian exponential function on a matrix group under a
%     left-invariant Riemannian metric.
%     The matrix group is defined by its Lie algebra, which is specified by
%     the columns of matBases. The metric is defined, on the algebra, by
%     the matrix MetricTensor. The starting point of the exponential is the
%     identity and U0 the initial velocity of the corresponding geodesic.
%     Q1 is the Riemannian exponential at t=1.
% 
%   [ Q1 , U1 ] = Exp( U0 , matBases , MetricTensor )
%     Returns in U1 the velocity in the algebra of the geodesic at t=1.
% 
%   [ Q1 , U1 , N_ITERS] = Exp( U0 , matBases , MetricTensor )
%     Returns in  N_ITERS(1)  the number of integration steps
%     and in  N_ITERS(2)  the total number of reductions.
% 
%   Qt = Exp( U0 , matBases , MetricTensor , TIMES )
%     TIMES is a vector of increasing and non-negative times where Qt is
%     computed.
% 
%   Qt = Exp( U0 , matBases , MetricTensor , OPTIONS )
%     OPTIONS is a structure with the following fields
%     (field names must be upper-case and if the field is not specified,
%     default parameters are used):
% 
%   OPTIONS.T           : TIMES where compute Qt and Ut (default  [ 1 ]).
%   OPTIONS.R           : Taylor order (default, depending on the initial
%                         velocity norm).
%                         ORDER is an scalar specifing the order R of the 
%                         Taylor expansion.
%                         If no ORDER is specified, the algorithm selects
%                         an order depending on the norm of the initial velocity.
%   OPTIONS.TOLERANCE   : maximum allowed norm of the last computed term in
%                         the Taylor expansion. It is used to calculate the
%                         value of delta to start the reducing process
%                         (default 1e-12).
%   OPTIONS.THRESHOLD_E : threshold to control the deviation of the
%                         conserved quantity E. If THRESHOLD_E is positive,
%                         step-size is reduced until 
%                                 abs( E - Ej ) < THRESHOLD_E
%                         If negative, step-size is reduced  until 
%                                 abs( E - Ej ) < THRESHOLD_E * abs( Ej )
%                         If 0 hard-rescaling of U to preserve E.
%                         If Inf ignore the control of E. 
%                         (default  -1e-12).
%   OPTIONS.THRESHOLD_S : threshold to control the deviation of the
%                         conserved quantity S. If positive, absolute
%                         deviation control. If negative, relative. 
%                         If Inf the control of S is ignored. 
%                         (default -1e-12)
%   OPTIONS.ALPHA       : reduction factor   (default 0.75)
%   OPTIONS.STEPSIZE    : If STEPSIZE is not equal to zero, the fixed step
%                         Taylor method is used. Therefore, no error
%                         controls are used. Neither step-size estimation
%                         nor step-size reduction are done. (default 0).
%   OPTIONS.PROJECT_U   : True or False, if true, after the U0-series
%                         integration, a re-projection on the algebra is
%                         done (default false)
%   OPTIONS.CHECK_Q     : A handle to a function of the form 
%                           @(Q) check_if_Q_belong_to_group(Q)
%                         and the function must return true or false. If
%                         this function is specified, after each Q-series
%                         integration step, Q is checked and if false a
%                         reduction of step-size is done. If empty, no 
%                         check is done. (default empty)
%   OPTIONS.PROJECT_Q   : A handle to a function of the form 
%                           Q_projected = @(Q) re_project_on_group(Q)
%                         If specified, after each Q-series integration
%                         step a reprojection is done. If empty, no
%                         re-projection is done. (default empty).
% 
%
%
%   Examples:
%     a) left-invariant geodesics on ST(1)
% 
%         matBases = [ 1 0 0 0 ; 0 0 1 0].';
%         MetricTensor = eye(2);
%         Q0 = [1 2;0 1];
%         V0 = [-2 5; 0 0];
%         Q1 = Exp( Q0 , V0 , matBases , MetricTensor );
% 
%     b) checking left-invariance
% 
%         Q0*Exp( [] , Q0 \ V0 , matBases , MetricTensor ) - Q1
% 
%     c) geodesic path
%
%         N = 1000;
%         Q = Q0; V = V0/N;
%         s = Q(1,1); t = Q(1,2);
%         for j = 1:N
%           [Q,V] = Exp( Q , V , matBases , MetricTensor );
%           s(j+1) = Q(1,1); t(j+1) = Q(1,2);
%         end
%         figure; plot( s , t ), axis equal
%
% 
%     d) checking function
% 
%         check_ST1 = @(Q) isequal( [ abs(Q(1,1)) , Q(1,2) ; 0 , 1 ] , Q )
%         Q1 = Exp( Q0 , V0 , matBases , MetricTensor , struct('CHECQ_Q',check_ST1) );
% 
%     e) PROJECT_Q function
% 
%         project_ST1 = @(Q) [ abs(Q(1,1)) , Q(1,2) ; 0 , 1 ];
%         Q1 = Exp( Q0 , V0 , matBases , MetricTensor , struct('PROJECT_Q',project_ST1) );
% 
%     f) changing MetricTensor by a scalar
% 
%         Exp( Q0 , V0 , matBases , MetricTensor ) - ...
%         Exp( Q0 , V0 , matBases , MetricTensor*100 )
% 
%        the result does not change if a scalar factor change MetricTensor
% 
%     g) changing MetricTensor by conjugacy
% 
%         C = [ 2 -3 ; 0 1 ];
%         MetricTensor_C = matBases.' * kron( inv(C) , C.' ) * pinv( matBases ).' * ...
%                          MetricTensor                                           * ...
%                          pinv( matBases ) * kron( inv( C ).' , C ) * matBases ;
%         Exp( Q0 , V0 , matBases , MetricTensor_C ) -...
%         inv( C ) * Exp( C*Q0*inv(C) , C*V0*inv(C) , matBases , MetricTensor ) * C
% 
%     h) computing Exp for a right-invariant Riemannian metric
% 
%         Exp_right = @(Q,V,mB,MT,varargin) inv( Exp( inv(Q) , -inv(Q)*V*inv(Q) , mB , MT ,varargin{:}) );
% 
%         Q1r = Exp_right( Q0 , V0 , matBases , MetricTensor )
% 
%     i) and check the right-invariance
% 
%         Exp_right( eye(2) , V0*inv(Q0) , matBases , MetricTensor )*Q0 - Q1r
% 
% 
% 
% Others bases to use:
% 
% matBases_ST1  = [1,0;zeros(1,2);0,1;zeros(1,2)];
% 
% matBases_ST2  = [1,0,0;zeros(1,3);zeros(1,3);zeros(1,3);1,0,0;zeros(1,3);0,1,0;0,0,1;zeros(1,3)];
% matBases_SO2  = [0;1;-1;0];
% matBases_SE2  = [zeros(1,3);1,0,0;zeros(1,3);-1,0,0;zeros(1,3);zeros(1,3);0,1,0;0,0,1;zeros(1,3)];
% matBases_SIM2 = [1,0,0,0;0,1,0,0;zeros(1,4);0,-1,0,0;1,0,0,0;zeros(1,4);0,0,1,0;0,0,0,1;zeros(1,4)];
% matBases_GL2  = eye([4,4]);
% matBases_SL2  = [1,0,0;0,1,0;0,0,1;-1,-0,-0];
% matBases_M2   = [zeros(1,6);0,1,0,0,0,0;0,0,1,0,-1,-0;0,0,-1,-0,-1,-0;0,-1,0,0,0,0;zeros(1,6);0,0,0,1,-0,-1;0,0,-0,-1,-0,-1;0,0,-1,-0,1,0;0,0,-0,-1,0,1;zeros(1,6);-1,0,0,0,0,0;0,0,-1,-0,-1,-0;0,0,-0,-1,-0,-1;-1,0,0,0,0,0;zeros(1,6)];
% matBases_GA2  = [1,0,0,0,0,0;0,0,0,1,0,0;zeros(1,6);0,1,0,0,0,0;0,0,0,0,1,0;zeros(1,6);0,0,1,0,0,0;0,0,0,0,0,1;zeros(1,6)];
% matBases_SA2  = [1,0,0,0,0;0,1,0,0,0;zeros(1,5);0,0,1,0,0;-1,-0,-0,0,0;zeros(1,5);0,0,0,1,0;0,0,0,0,1;zeros(1,5)];
% 
% matBases_ST3  = [1,0,0,0;zeros(1,4);zeros(1,4);zeros(1,4);zeros(1,4);1,0,0,0;zeros(1,4);zeros(1,4);zeros(1,4);zeros(1,4);1,0,0,0;zeros(1,4);0,1,0,0;0,0,1,0;0,0,0,1;zeros(1,4)];
% matBases_SO3  = [zeros(1,3);1,0,0;0,1,0;-1,0,0;zeros(1,3);0,0,1;0,-1,0;0,0,-1;zeros(1,3)];
% matBases_SE3  = [zeros(1,6);1,0,0,0,0,0;0,1,0,0,0,0;zeros(1,6);-1,0,0,0,0,0;zeros(1,6);0,0,1,0,0,0;zeros(1,6);0,-1,0,0,0,0;0,0,-1,0,0,0;zeros(1,6);zeros(1,6);0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;zeros(1,6)];
% matBases_SIM3 = [1,0,0,0,0,0,0;0,1,0,0,0,0,0;0,0,1,0,0,0,0;zeros(1,7);0,-1,0,0,0,0,0;1,0,0,0,0,0,0;0,0,0,1,0,0,0;zeros(1,7);0,0,-1,0,0,0,0;0,0,0,-1,0,0,0;1,0,0,0,0,0,0;zeros(1,7);0,0,0,0,1,0,0;0,0,0,0,0,1,0;0,0,0,0,0,0,1;zeros(1,7)];
% matBases_GL3  = eye([9,9]);
% matBases_SL3  = [1,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,1,0,0,0;0,0,0,0,0,1,0,0;zeros(1,6),1,0;zeros(1,7),1;-1,-0,-0,-0,-1,-0,-0,-0];
% matBases_M3   = [zeros(1,10);0,1,zeros(1,8);0,0,1,0,0,0,0,0,0,0;0,0,0,0,1,0,0,-1,-0,-0;0,0,0,0,-1,-0,-0,-1,-0,-0;0,-1,zeros(1,8);zeros(1,10);0,0,0,1,0,0,0,0,0,0;0,0,0,0,0,1,0,-0,-1,-0;0,0,0,0,-0,-1,-0,-0,-1,-0;0,0,-1,0,0,0,0,0,0,0;0,0,0,-1,0,0,0,0,0,0;zeros(1,10);zeros(1,6),1,-0,-0,-1;0,0,0,0,-0,-0,-1,-0,-0,-1;0,0,0,0,-1,-0,-0,1,0,0;0,0,0,0,-0,-1,-0,0,1,0;0,0,0,0,-0,-0,-1,0,0,1;zeros(1,10);-1,zeros(1,9);0,0,0,0,-1,-0,-0,-1,-0,-0;0,0,0,0,-0,-1,-0,-0,-1,-0;0,0,0,0,-0,-0,-1,-0,-0,-1;-1,zeros(1,9);zeros(1,10)];
% matBases_GA3  = [1,zeros(1,11);0,0,0,0,1,0,0,0,0,0,0,0;zeros(1,8),1,0,0,0;zeros(1,12);0,1,zeros(1,10);0,0,0,0,0,1,0,0,0,0,0,0;zeros(1,9),1,0,0;zeros(1,12);0,0,1,zeros(1,9);zeros(1,6),1,0,0,0,0,0;zeros(1,10),1,0;zeros(1,12);0,0,0,1,zeros(1,8);zeros(1,7),1,0,0,0,0;zeros(1,11),1;zeros(1,12)];
% matBases_SA3  = [1,zeros(1,10);0,1,zeros(1,9);0,0,1,zeros(1,8);zeros(1,11);0,0,0,1,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0,0;zeros(1,11);zeros(1,6),1,0,0,0,0;zeros(1,7),1,0,0,0;-1,-0,-0,-0,-1,-0,-0,-0,0,0,0;zeros(1,11);zeros(1,8),1,0,0;zeros(1,9),1,0;zeros(1,10),1;zeros(1,11)];
% 

%   Reference:
%   Zacur E., Bossa M., Olmos S., "A  method to compute left-invariant
%   Riemannian geodesics on matrix groups. Illustration on spatial 
%   transformations"
%
%   Ernesto Zacur
%   Universidad de Zaragoza
%   $Date: 2012/10/09$


  if nargin < 3
    error('Exp:NoEnoughArguments','At least 3 arguments were expected.');
  end

  %defaults
  R           = [];       %if empty, use a rule of thumb to decide the order
  T           = 1;        %integration time
  TOLERANCE   =  1e-10;
  THRESHOLD_E = -1e-10;
  THRESHOLD_S = -1e-10;
  ALPHA       = 0.75;
  STEPSIZE    = 0;
  MINSTEPSIZE = 1e-4;
  PROJECT_U   = false;
  CHECK_Q     = [];
  PROJECT_Q   = [];
  VERBOSE     = false;
  
  if nargin > 3 
    if isstruct( OPTS )
      if isfield( OPTS , 'T'           ), T           = OPTS.T;           OPTS = rmfield(OPTS,'T');            end
      if isfield( OPTS , 'ORDER'       ), R           = OPTS.ORDER;       OPTS = rmfield(OPTS,'ORDER');        end
      if isfield( OPTS , 'TOLERANCE'   ), TOLERANCE   = OPTS.TOLERANCE;   OPTS = rmfield(OPTS,'TOLERANCE');    end
      if isfield( OPTS , 'THRESHOLD_E' ), THRESHOLD_E = OPTS.THRESHOLD_E; OPTS = rmfield(OPTS,'THRESHOLD_E');  end
      if isfield( OPTS , 'THRESHOLD_S' ), THRESHOLD_S = OPTS.THRESHOLD_S; OPTS = rmfield(OPTS,'THRESHOLD_S');  end
      if isfield( OPTS , 'ALPHA'       ), ALPHA       = OPTS.ALPHA;       OPTS = rmfield(OPTS,'ALPHA');        end
      if isfield( OPTS , 'STEPSIZE'    ), STEPSIZE    = OPTS.STEPSIZE;    OPTS = rmfield(OPTS,'STEPSIZE');     end
      if isfield( OPTS , 'MINSTEPSIZE' ), MINSTEPSIZE = OPTS.MINSTEPSIZE; OPTS = rmfield(OPTS,'MINSTEPSIZE');  end
      if isfield( OPTS , 'PROJECT_U'   ), PROJECT_U   = OPTS.PROJECT_U;   OPTS = rmfield(OPTS,'PROJECT_U');    end
      if isfield( OPTS , 'CHECK_Q'     ), CHECK_Q     = OPTS.CHECK_Q;     OPTS = rmfield(OPTS,'CHECK_Q');      end
      if isfield( OPTS , 'PROJECT_Q'   ), PROJECT_Q   = OPTS.PROJECT_Q;   OPTS = rmfield(OPTS,'PROJECT_Q');    end
      if isfield( OPTS , 'VERBOSE'     ), VERBOSE     = OPTS.VERBOSE;     OPTS = rmfield(OPTS,'VERBOSE');      end
      if ~isempty( fieldnames(OPTS) )
        warning('Exp:unknowOPTS', 'Unknown OPTIONS in OPTS. Check all OPTIONS are specified as upper-case.');
      end
      
    elseif isnumeric( OPTS )
      T = OPTS;
    else
      error('Exp:IncorrectArguments','A struct or a vector was expected.');
    end
  end

  if ~isvector( T ) || min( T ) < 0 ||  any( diff(T) <= 0 )
    error('Exp:InvalidIntegrationTime','T is expected to be an increasing and non-negative vector.');
  end
  T = T(:).';
  
  if ~isscalar( STEPSIZE ) || STEPSIZE < 0
    error('Exp:BadSTEPSIZE', 'STEPSIZE must be a non-negative scalar.');
  end
  if STEPSIZE ~= 0       %in case of stepsize ~= 0, fixed step algorithm will be switched
                         %therefore, no reductions will be done
    THRESHOLD_E = Inf;
    THRESHOLD_S = Inf;
    CHECK_Q     = [];
  end

  if ~isscalar( THRESHOLD_E )
    error('Exp:BadTHRESHOLD_E', 'THRESHOLD_E must be a scalar.');
  end
  if isfinite(THRESHOLD_E) && THRESHOLD_E < 0 && abs( THRESHOLD_E ) < 100 * eps(1)
    THRESHOLD_E = -100 * eps(1);
    warning('Exp:RelThrEIncrease', 'THRESHOLD_E has been increased to -100*eps(1)')
  end
  if ~isscalar( THRESHOLD_S ) || isequal( THRESHOLD_S , 0 )
    error('Exp:BadTHRESHOLD_S', 'THRESHOLD_S must be a nonzero scalar.');
  end
  if isfinite(THRESHOLD_S) && THRESHOLD_S < 0 && abs( THRESHOLD_S ) < 100 * eps(1)
    THRESHOLD_S = -100 * eps(1);
    warning('Exp:RelThrSIncrease', 'THRESHOLD_S has been increased to -100*eps(1)')
  end
  if MINSTEPSIZE < 0
    error('Exp:BadMINSTEPSIZE','MINSTEPSIZE can not be negative');
  end
  if MINSTEPSIZE < 100*eps(1)
    MINSTEPSIZE =  100 * eps(1);
    warning('Exp:MINSTEPSIZEIncrease', 'MINSTEPSIZE has been increased to 100*eps(1)')
  end
  
  %checking matrices dimension
  n  = size( U0  , 1 );
  k  = size( B  , 2 );
  if ~isequal( size(U0)             , [n   , n]), error('Exp:IncorrectU0','A square matrix was expected in U0'); end
  if ~isequal( size(B )             , [n*n , k]), error('Exp:IncorrectB','Basis matrix have to be  (n*n) x k'); end
  if ~isequal( size(MetricTensor )  , [k   , k]), error('Exp:IncorrectMetricTensor','Metric tensor have to be  k x k'); end
  Q0 = eye(n,n); Qtt = eye(n);
  
  if ~isempty( CHECK_Q )
    try
      checkQ = feval( CHECK_Q , Q0 );
    catch LE
      if ~isempty(which('disperror')),   disperror( LE ); end
      error('Exp:errorInCHECK_Q','Invalid CHECK_Q function.');
    end
    if ~checkQ
      error('Exp:I_NotBelongToGroup','Identity belong to the group. It must satisfy CHECK_Q');
    end
  end
  if PROJECT_U
    PROJECTION_MATRIX = B * ( B \ eye(n*n,n*n) );
  end
  
  
  %linear operator sigma
  SIGMA = B * ( MetricTensor \ B.' );
  if any( isnan( SIGMA(:) ) )
    SIGMA = B * pinv( MetricTensor ) * B.';
  end

  %linear operator chi
%   Z = null( B.' );
%   CHI   = ( SIGMA + Z*Z.' ) \ eye(n*n,n*n);
  CHI = B \ eye(n*n,n*n);
  CHI = CHI.' * MetricTensor * CHI;

  %some usefulls
  vec  = @(x) x(:);
  ITERATIONS = [0 0];    

  %outputs
  Qt = zeros(n,n,numel(T));
  Ut = zeros(n,n,numel(T));
  lastT = 0;
  
  
  %check uniparametric condition
  xU0 = reshape( CHI*U0(:) ,[n,n] );
  if VERBOSE
    fprintf('to compute as UNIPARAMETRIC: %g\n', max( abs( SIGMA * vec( U0.'*xU0 - xU0 * U0.' ) ) ) );
  end
  if  STEPSIZE == 0  &&  max( abs( SIGMA * vec( U0.'*xU0 - xU0 * U0.' ) ) ) < 1e-14
    %checking if it fulfill  SIGMA[ U0.' , CHI( U0 ) ] == 0
    
    if     isequal( U0 , diag(diag(U0)) )                                     %scales case
      for tt = 1:numel(T),  Qt(:,:,tt) = diag( exp( diag(U0) * T(tt) ) );  end
    elseif isequal( U0 , [ zeros(n-1,n-1) , ones(n-1,1) ; zeros(1,n) ].*U0 )  %translations case
      for tt = 1:numel(T),  Qt(:,:,tt) = eye(n) + U0*T(tt);                end
    else                                                                      %general case
      for tt = 1:numel(T),  Qt(:,:,tt) = expm( U0*T(tt) );                 end
    end
    
    Ut = bsxfun(@plus,U0,Ut);
    
    if VERBOSE, fprintf('COMPUTED AS UNIPARAMETRIC!\n'); end
    
  else
  
%     if isfinite( THRESHOLD_S )
%       condCHI = cond(CHI);
%       if VERBOSE, fprintf('cond( CHI ) = %g\n',condCHI ); end
%       if condCHI > 1e2
%         warning('Exp:BadConditionedCHI', 'The operator CHI is bad conditioned.\nPossible inacurracies in step reduction process.');
%       end
%     end

    ZERO = zeros(n,n);
    
    
    %rule of thumb to determine R
    if isempty(R)
      mu = B \ U0(:);
      normU0 = sqrt( mu.' * MetricTensor * mu );
      R = ceil( 15*realpow( normU0 , 1/6 ) );
    end

    if VERBOSE, fprintf('computation using R: %d\n\n', R ); end

    %precompute BINOMIAL coefficients
    BINOMIAL   = [1,zeros(1,11);1,1,zeros(1,10);1,2,1,zeros(1,9);1,3,3,1,zeros(1,8);1,4,6,4,1,0,0,0,0,0,0,0;1,5,10,10,5,1,0,0,0,0,0,0;1,6,15,20,15,6,1,0,0,0,0,0;1,7,21,35,35,21,7,1,0,0,0,0;1,8,28,56,70,56,28,8,1,0,0,0;1,9,36,84,126,126,84,36,9,1,0,0;1,10,45,120,210,252,210,120,45,10,1,0;1,11,55,165,330,462,462,330,165,55,11,1];
    if R+1 > size(BINOMIAL,1), BINOMIAL(R+1,R+1) = 0; end
    for j = 12:R
      BINOMIAL(j+1,:) = [ 1 , BINOMIAL( j , 1:end-1 ) + BINOMIAL( j , 2:end ) ];
    end

    %precompute usefull constants to calculate delta*
    R_factorial = prod( 1:R );
    inv_R       = 1/R;


    %containers to store U0.(m), chi( U0.(m) ) and tau( U0.(m) )
    Us = cell(R+1,1);
    xU = cell(R  ,1);
    tU = cell(R  ,1);
    %container to store Q.(m)
    Qs = cell(R+1,1);
    Qs{1} = eye(n,n);


    t = 0;
    if VERBOSE, fprintf('\n*****************\nt: %5.20g\n', t ); end
    
    if T(1) == t
      Qt(:,:,1) = Q0*Qtt;
      Ut(:,:,1) = U0;
    end
    
    while t < T(end)
      E0 = trace( xU0.' * U0 );
      S0 = U0(:);
      
      Us{1} = U0;
      xU{1} = xU0;
      tU{1} = U0.';
      
      %Qs{1} = eye(n);

      %%recursive process to compute the taylor coefficients
      for m = 2:R
        Us{m} = ZERO;
        for j=1:m-1
          Us{m} = Us{m} + BINOMIAL(m-1,j)* reshape( SIGMA * vec( tU{j}*xU{m-j}  -  xU{m-j}*tU{j}  ) ,[n,n]);
        end
        %se podria proyectar aca U para no propagar los errores.?
        xU{m} = reshape( CHI*Us{m}(:) ,[n,n]);
        tU{m} =    Us{m}.';

        Qs{m} = ZERO;
        for j=1:m-1
          Qs{m} = Qs{m} + BINOMIAL(m-1,j) * Qs{j}*Us{m-j};
        end
      end

      Us{R+1} = ZERO;
      for j=1:R
        Us{R+1} = Us{R+1} + BINOMIAL(R,j)* reshape( SIGMA * vec( tU{j}*xU{R+1-j}  -  xU{R+1-j}*tU{j} ) ,[n,n]);
      end

      Qs{R+1} = ZERO;
      for j=1:R
        Qs{R+1} = Qs{R+1} + BINOMIAL(R,j)* Qs{j}*Us{R+1-j};
      end
      %%END recursive process to compute the taylor coefficients

      
      if STEPSIZE == 0
        
        %estimate delta
        delta = min( T(end)-t  ,...
                     realpow( R_factorial * TOLERANCE / max( max(abs( Qs{R+1}(:) )) ,...
                                                             max(abs( Us{R+1}(:) )) ), inv_R ) );
        if VERBOSE, fprintf('delta* to start reduction: %5.20g\n', delta ); end
        EPSt = max( 100*eps( t ) , MINSTEPSIZE );

      else
        
        delta = min( STEPSIZE , T(end)-t );
        if VERBOSE, fprintf('fixed step delta: %5.20g\n', delta ); end

      end
      
      %reduction process if needed
      while 1
        [U0,Qd] = taylorSerie( delta );
        
        %project onto the algebra
        %posiblemente lo mejor seria proyectar despues de hacer las
        %reducciones necesarias, es que asi, da peor!!!
        if PROJECT_U
          Up = reshape( PROJECTION_MATRIX * U0(:) ,[n,n]);
          if VERBOSE, fprintf('difference with projected U0: %5.20g\n', max( abs( U0(:) - Up(:) ) ) ); end
          U0 = Up;
        end

        %rescaled on the ellipse
        %posiblemente lo mejor seria rescalar despues de hacer las
        %reducciones necesarias, es que asi, da peor!!!
        if THRESHOLD_E == 0
          Ur = U0;
          f  = E0/trace( reshape( CHI*Ur(:) ,[n,n]).' * Ur  );
          for rescalings = 1:5
            if f == 1, break; end
            Ur = Ur*f;
            f  = E0/trace( reshape( CHI*Ur(:) ,[n,n]).' * Ur  );
          end
          Ur = Ur*f;
          if VERBOSE, fprintf('difference with rescaled U0: %5.20g\n', max( abs( U0(:) - Ur(:) ) ) ); end
          U0 = Ur;
        end

        %project onto the group
        if ~isempty( PROJECT_Q )
          try
            Qp = feval( PROJECT_Q , Qd );
          catch LE
            if ~isempty(which('disperror')),   disperror( LE ); end
            error('Exp:errorInPROJECT_Q','Error in PROJECT_Q function.');
          end
          if VERBOSE, fprintf('difference with projected Qd: %5.20g\n', max( abs( Qd(:) - Qp(:) ) ) ); end
          Qd = Qp;
        end

        xU0 = reshape( CHI*U0(:) ,[n,n]);

        if STEPSIZE ~= 0, break; end    %no reductions in case of fixed step
        
        checkQ = true;
        if ~isempty( CHECK_Q )
          try
            checkQ = feval( CHECK_Q , Qd );
          catch LE
            if ~isempty(which('disperror')),   disperror( LE ); end
            error('Exp:errorInCHECK_Q','Error in CHECK_Q function.');
          end
        end
        
        if ~checkQ
          if VERBOSE, fprintf( 'Q3 does not satisfy CHECK_Q\n' ); end
        else
          %if THRESHOLD is greater than zero, consider it as absolute
          %threshold, else, conside it as relative threshold
          %if it is inf, ignore

          if  isfinite( THRESHOLD_E ), E = trace( xU0.' * U0 ); end
          if ~isfinite( THRESHOLD_E )                                        ||...
               THRESHOLD_E == 0                                              ||...
             ( THRESHOLD_E > 0 && abs( E - E0 ) <=      THRESHOLD_E        ) ||...
             ( THRESHOLD_E < 0 && abs( E - E0 ) <= abs( THRESHOLD_E * E0 ) ) 
               
            if  isfinite( THRESHOLD_S ), S = reshape( SIGMA * vec( ( Qd.' ) \ ( xU0 * Qd.' ) ) ,[n,n]); end
            if ~isfinite( THRESHOLD_S )                                                  ||...
               ( THRESHOLD_S > 0 && all( abs( S(:) - S0 ) <=      THRESHOLD_S        ) ) ||...
               ( THRESHOLD_S < 0 && all( abs( S(:) - S0 ) <= abs( THRESHOLD_S * S0 ) ) )

             %control conditions are satisfied, reduction process already done
              break

            elseif VERBOSE, fprintf( '| S - Sj |: %5.20g\n', max( abs( S(:) - S0 ) ) );
            end
          elseif VERBOSE, fprintf( '| E - Ej |: %5.20g\n', max( abs( E - E0 ) ) );
          end
        end

        if delta < EPSt
          warning('Exp:TooSmallStep', 'The reduction process reachs to a very small step.\nUnable to meet integration thresholds ( at t=%1.15g ).\nContinue without controls in this step.', t );
          break;
        end

        %reduce STEP-SIZE
        delta = delta * ALPHA;
        ITERATIONS(2) = ITERATIONS(2) + 1;
        if VERBOSE, fprintf('reducing step to: %5.20g\n', delta ); end
      end


      %fill outputs
      for tt = find( T > t & T <= t+delta )
        [Utt,Qtt] = taylorSerie( T(tt) - t );
        Qt(:,:,tt) = Q0*Qtt;
        Ut(:,:,tt) = Utt;
      end
      %END fill outputs
      
      
      %update step
      ITERATIONS(1) = ITERATIONS(1) + 1;
      t = t + delta;
      Q0 = Q0*Qd;
      if VERBOSE, fprintf('\n*****************\nt: %5.20g\n', t ); end
    end
    
  end

  if VERBOSE, fprintf('\n-----------------\niterations: %d\nreductions: %d\n\n', ITERATIONS ); end

  
  function [U,Q] = taylorSerie( delta )
    f = 1;
    U = Us{1};
    Q = Qs{1};
    for mm = 2:R+1
      f = f * delta/(mm-1);
      U = U + f*Us{mm};
      Q = Q + f*Qs{mm};
    end
  end

end
