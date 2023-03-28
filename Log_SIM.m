function [d,U] = Log_SIM( Q )

  d = size( Q , 1 );
  d = d-1;

  S = det( Q(1:d,1:d) )^(1/d);              R = Q(1:d,1:d)/S;         T = Q(1:d,d+1);
  if true
    if ~isequal( size(Q) , [d+1 d+1] ), error('Q must be square'); end
    
    if any( Q(d+1,1:d) ),   error('la ultima fila no es zero'); end
    if Q(d+1,d+1) ~= 1,     error('el ultimo elemento no es 1'); end
    
    Q_ = [ S*R , T ; zeros(1,d) , 1 ];
    if maxnorm( Q , Q_ ) > 1e-10,             error('bad structured Q');        end
    
    if maxnorm( R   * R.' - eye(d) ) > 1e-10, error('no satisface R*R'' == I'); end
    if maxnorm( R.' * R   - eye(d) ) > 1e-10, error('no satisface R''*R == I'); end
  end


  ST = [ eye(d)*S , T ; zeros(1,d) , 1 ];

  if nargout < 2
  
    dSO = Log_SO( R  );
    dST = Log_ST( ST );
    
  else    

    [dSO,VSO] = Log_SO( R  );
    [dST,VST] = Log_ST( ST );
    
    U = zeros(d+1,d+1);
      
    vS0 = VSO(:,:);
    vS0(d+1,d+1) = 0;
    vST = VST(:,:);

    U( : , : ) = ( vS0 + vST );
    
  end

  
  d = sqrt( dSO.^2 + dST.^2 );
  
end
