function [B,m] = LIEbasis( B , m )

%{
  
B = LIEbasis('sim(2)');
K = size(B,2);
m = eye(K)+2*randn(K); m = m*m.';

B1 = LIEbasis(B,m);
maxnorm( pinv(B).'*m*pinv(B)  - pinv(B1).'*pinv(B1) )

[Bo,mm]= LIEbasis( B , m );
maxnorm( pinv(B).'*m*pinv(B)  - pinv(Bo).'* mm * pinv(Bo) )
maxnorm( eye(K) - Bo.'*Bo )

[Bo,mm]= LIEbasis( LIEbasis(B,m) );
maxnorm( pinv(B).'*m*pinv(B)  - pinv(Bo).'* mm * pinv(Bo) )
maxnorm( eye(K) - Bo.'*Bo )
  
%}



  CHECK_BASIS = true;
  if ischar( B )
    CHECK_BASIS = false;
    B = lower( B );
    if 0
    elseif ~isempty( regexp( B , '^t\((\d+)\)' , 'once' ) )
      N  = regexp( B , 't\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 't(n)';
    elseif ~isempty( regexp( B , '^s\((\d+)\)' , 'once' ) )
      N  = regexp( B , 's\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 's(n)';
    elseif ~isempty( regexp( B , '^s\+t\((\d+)\)' , 'once' ) )
      N  = regexp( B , 's\+t\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 's+t(n)';
    elseif ~isempty( regexp( B , '^st\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'st\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'st(n)';
    elseif ~isempty( regexp( B , '^so\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'so\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'so(n)';
    elseif ~isempty( regexp( B , '^so\*s\((\d+)\)' , 'once' ) )
      N  = regexp( B , '^so\*s\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'so*s(n)';
    elseif ~isempty( regexp( B , '^so\*t\((\d+)\)' , 'once' ) )
      N  = regexp( B , '^so\*t\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'so*t(n)';
    elseif ~isempty( regexp( B , '^so\*s\*t\((\d+)\)' , 'once' ) )
      N  = regexp( B , '^so\*s\*t\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'so*s*t(n)';
    elseif ~isempty( regexp( B , '^s\*so\*t\((\d+)\)' , 'once' ) )
      N  = regexp( B , '^s\*so\*t\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 's*so*t(n)';
    elseif ~isempty( regexp( B , '^se\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'se\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'se(n)';
    elseif ~isempty( regexp( B , '^sim\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'sim\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'sim(n)';
    elseif ~isempty( regexp( B , '^sl\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'sl\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'sl(n)';
    elseif ~isempty( regexp( B , '^sa\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'sa\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'sa(n)';
    elseif ~isempty( regexp( B , '^gl\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'gl\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'gl(n)';
    elseif ~isempty( regexp( B , '^ga\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'ga\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'ga(n)';
    elseif ~isempty( regexp( B , '^mo\((\d+)\)' , 'once' ) )
      N  = regexp( B , 'mo\((\d+)\)' ,'tokens' );
      N  = str2double( N{1} );
      B = 'mo(n)';
    end

    
    switch B
      case 'mo(n)'
%         if N ~= 2, error('se esperaba que n sea 2'); end
%         B = { ...
%               kron( [1 0;0 -1] , [1 0; 0 1] ) ,...
%               kron( [0 1;0  0] , [1 0; 0 1] ) ,...
%               kron( [0 0;1  0] , [1 0; 0 1] ) ,...
%               kron( [1 0;0 -1] , [0 1;-1 0] ) ,...
%               kron( [0 1;0  0] , [0 1;-1 0] ) ,...
%               kron( [0 0;1  0] , [0 1;-1 0] ) };

        B = {};
        
        B{end+1} = blkdiag( zeros(N) , [0 -1; -1 0] );

        bi = zeros( N*(N-1)/2 , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = blkdiag( skewmatrix( bi ) , zeros(2,2) );
          bi(i) = 0;
        end

        bi = zeros( N , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = [ [ zeros(N,N) , -bi , -bi ] ; [ [ bi.' ; -bi.' ] , zeros(2,2) ] ];
          bi(i) = 0;
        end
        
        bi = zeros( N , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = [ [ zeros(N,N) , bi , -bi ] ; [ [ -bi.' ; -bi.' ] , zeros(2,2) ] ];
          bi(i) = 0;
        end
      
      case 't(n)'
        B = {};
        
        bi = zeros(N+1,N+1);
        for i = 1:N
          bi(i,N+1) = 1;
          B{end+1} = bi;
          bi(i,N+1) = 0;
        end

      case 's(n)'
        B = {};
        
        bi = zeros(N,1);
        for i = 1:N
          bi(i) = 1;
          B{end+1} = diag(bi);
          bi(i) = 0;
        end

      case 's+t(n)'
        B = {};
        
        bi = zeros( N+1 , N+1 );
        for i = 1:N
          bi(i,i)  = 1;
          B{end+1} = bi;
          bi(i,i)  = 0;
        end
        
        for i = 1:N
          bi(i,N+1) = 1;
          B{end+1} = bi;
          bi(i,N+1) = 0;
        end
%         B = cellfun( @(b) b/norm(b,'fro') , B , 'UniformOutput',false );

      case 'st(n)'
        B = {};
        
        bi = eye( N+1 , N+1 );
        bi(N+1,N+1) = 0;
        B{end+1} = bi;
        
        bi = zeros( N+1 , N+1 );
        for i = 1:N
          bi(i,N+1) = 1;
          B{end+1}  = bi;
          bi(i,N+1) = 0;
        end
%         B = cellfun( @(b) b/norm(b,'fro') , B , 'UniformOutput',false );

      case 'so(n)'
        B = {};

        bi = zeros( N*(N-1)/2 , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = skewmatrix( bi );
          bi(i) = 0;
        end
        
      case 'so*s(n)'
        B = {};

        bi = zeros( N*(N-1)/2 , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = skewmatrix( bi );
          bi(i) = 0;
        end

        B{end+1} = eye(N);

      case 'so*t(n)'
        B = {};

        bi = zeros( N*(N-1)/2 , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = blkdiag( skewmatrix( bi ) , zeros(N+1) );
          bi(i) = 0;
        end

        bi = zeros( N , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = blkdiag( zeros(N) , [ zeros(N,N) , bi ; zeros(1,N) , 0 ] );
          bi(i) = 0;
        end
        
      case 'so*s*t(n)'
        B = {};

        bi = zeros( N*(N-1)/2 , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = blkdiag( skewmatrix( bi ) , zeros(N+1) );
          bi(i) = 0;
        end

        B{end+1} = blkdiag( eye(N) , zeros(N+1) );
        
        bi = zeros( N , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = blkdiag( zeros(N) , [ zeros(N,N) , bi ; zeros(1,N) , 0 ] );
          bi(i) = 0;
        end
        
      case 's*so*t(n)'
        B = {};

        bi = zeros( N*(N-1)/2 , 1 );
        B{end+1} = blkdiag( eye(N) , zeros(N+1) );

        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = blkdiag( skewmatrix( bi ) , zeros(N+1) );
          bi(i) = 0;
        end

        
        bi = zeros( N , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = blkdiag( zeros(N) , [ zeros(N,N) , bi ; zeros(1,N) , 0 ] );
          bi(i) = 0;
        end
        
      case 'se(n)'
        B = {};

        bi = zeros( N*(N-1)/2 , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = [ skewmatrix( bi ) , zeros(N,1) ; zeros(1,N) , 0 ];
          bi(i) = 0;
        end

        bi = zeros( N , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = [ zeros(N,N) , bi ; zeros(1,N) , 0 ];
          bi(i) = 0;
        end

      case 'sim(n)'
        B = {};
        
        s = eye( N+1 , N+1 );
        s(end) = 0;
        B{end+1} = s;

        bi = zeros( N*(N-1)/2 , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = [ skewmatrix( bi ) , zeros(N,1) ; zeros(1,N) , 0 ];
          bi(i) = 0;
        end

        bi = zeros( N , 1 );
        for i = 1:numel(bi);
          bi(i) = 1;
          B{end+1} = [ zeros(N,N) , bi ; zeros(1,N) , 0 ];
          bi(i) = 0;
        end

      case 'gl(n)'
        B = eye(N*N);

      case 'sl(n)'
        B = eye(N*N);
        
        bi = eye(N,N);
        bi(end,end) = 0;
        B(end,:) = - bi(:).';
        B(:,end) = [];
        
      case 'ga(n)'
        B = {};
        
        bi = zeros(N+1,N+1);
        for i = 1:N
          for j = 1:N+1
            bi(i,j)  = 1;
            B{end+1} = bi;
            bi(i,j)  = 0;
          end
        end

      case 'sa(n)'
        B = eye(N*N);
        
        bi = eye(N,N);
        bi(end,end) = 0;
        B(end,:) = - bi(:).';
        B(:,end) = [];
        
        B = arrayfun( @(j)tosquare(B(:,j)) , 1:size(B,2) , 'UniformOutput', false );
        add_0 = @(x) [ x  , zeros(N,1) ; zeros(1,N) , 0 ];
        B = cellfun( @(b) add_0(b) , B , 'UniformOutput', false );
        
        bi = zeros(N+1,N+1);
        for i = 1:N
          bi(i,N+1) = 1;
          B{end+1} = bi;
          bi(i,N+1) = 0;
        end
        
      otherwise
        error('unknow basis');  
    end

  end
  
  
  if iscell( B )
    if CHECK_BASIS
      issquare = @(x) (ndims(x)==2) && (size(x,1)==size(x,2));
      if any( ~cellfun( @(x) issquare(x) , B ) )
        error('the basis must be squares');
      end
    end
    B = cell2mat( cellfun( @(b) b(:) , B(:).' , 'UniformOutput', false ) );
  end
  
  
  N = sqrt( size(B,1) );
  K = size(B,2);
  
  
  if CHECK_BASIS
    if mod(N,1), error('no valid basis (no square size)'); end

    %check the rank
    if rank(B) ~= K, error('not all basis independent'); end

    %check if B is a Lie-Algebra
    for i = 1:K
      for j = i+1:K
        Bi = reshape( B(:,i) , [N,N] );
        Bj = reshape( B(:,j) , [N,N] );
        COMM = Bi*Bj - Bj*Bi;
        if rank( [ B , COMM(:) ] ) ~= K
          error('basis  %d  and  %d  generate a new independent algebra',i,j);
        end
      end
    end
  end

  
  if nargout > 1
    
    if nargin < 2 || isempty( m )
      m = eye(K);
    end
    
    [q,r] = qr( B );
    
    B = q(:,1:K);
    r = r(1:K,1:K);
    
    r = r \ eye(K);
    
    m = r.' * m * r;
    
  elseif ~( nargin < 2 || isempty( m ) || isequal( m , eye(K) ) )
    
    if ~isequal( size(m) , [K K] )
      error('incorrect metric tensor');
    end
%     try
%       chol( m );
%     catch
%       error('metric tensor is expected to be positive definite');
%     end
      
    [W,D] = eig( m );
    D = diag( D ).';
    if any( D <= 0 ) || any( imag(D) > 0 )
      error('metric tensor is expected to be symmetric and positive definite');
    end
    
    W = bsxfun( @rdivide , W , sqrt( D ) );
    B = B * W;
    
  end
  
  

end
