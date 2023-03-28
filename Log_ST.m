function [dist,U] = Log_ST( Q )

  d = size( Q , 1 );
  d = d-1;


  S = Q(1,1);    T = Q(1:d,d+1);

  if true
    if S <= 0
      error('S0 must be greater than zero.'); end

    Q_ = [ S*eye(d) , T ; zeros(1,d) , 1 ];

    if ~isequal( Q , Q_ ), error('bad structured X0'); end
  end
  
  Tn = sqrt( sum( T.^2 , 1 ) );
  

  dist = acosh( 1 + ( Tn.^2 + (S-1)^2 )./( 2*S ) );

  if nargout > 1
    
      if Tn ~= 0
        c = ( Tn^2 + S^2 - 1 )/( 2 * Tn );
        v = [ c ; T/Tn ];
        v = v/norm(v)*dist;
        
        
        U = blkdiag( eye(d)*v(1) , 0 );
        U(1:d,d+1) = v(2:end);
      else
        U = blkdiag( eye(d)*log( S ) , 0 );
      end

  end
  
end
