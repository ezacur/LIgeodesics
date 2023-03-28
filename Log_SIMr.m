function [d,U] = Log_SIMr( Q )

  [d,U] = Log_SIM( Q\eye(size(Q,1)) );
  U = -U;
  
end
