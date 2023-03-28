function Qt = Exp_ST1( U0 , time , M )

  if nargin < 1, error('at least 1 argument is expected'); end
  if ~isequal( size(U0) , [2 2] )
    error('U0 must be 2x2');
  end
  if nargin < 3
    M = [];
  end
  if nargin == 2 && isequal(  size( time ) , [2 2] )
    M = time;
    time = [];
  end
  if nargin < 2
    time = [];
  end
  if isempty( M    ), M    = eye(2); end
  if isempty( time ), time = 1;      end
  if ~isvector( time )
    error('time is expected to be a vector');
  end
  if ~isequal( M , 'g' )  &&  ~isequal( size(M) , [2 2] )
    error('M must be 2x2');
  end
  if max(abs( vec( M - M.' ) ) ) > 1e-10
    error('M must be symmetric');
  end


  if isequal( M , 'g' )
    Qt = zeros( 2,2, numel(time) );
    for k = 1:size(Qt,3)
      Qt(:,:,k) = expm( U0 * time(k) );
    end
    
    return;
  end
  
  

  s =  sqrt( det(M) )/M(2,2);
  t =  M(1,2)/M(2,2);
  
  Y = [s , t;0 , 1]; iY = [1/s , -t/s ; 0 , 1];
  
  U0 = iY * U0 * Y;
  
  Qt = Exp_ST( U0 , time );
  
  for k = 1:size(Qt,3)
    Qt(:,:,k) = Y * Qt(:,:,k) * iY;
  end
  
  
end
