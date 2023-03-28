function M = projectOnSL( M )

  sign_det = det(M);
  if sign_det == 1, return; end
  
  sign_det = sign( sign_det );
  [U,s,V] = svd( M );
  
  s = diag(s); s = s(:).';
  
%   norm2 = @(x) x(:).'*x(:);
%   s_opt = Optimize( @(ss) norm2( [ ss , sign_det/prod(ss) ] - s ) , s(1:end-1) ,...
%     'methods',{'conjugate','coordinate',1},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'i'}}),'verbose',0,'noplot');
  s_opt = Optimize( @(ss) ener(ss) , s(1:end-1) ,...
    'methods',{'conjugate','coordinate',1},'verbose',0,'noplot');

  M = U * diag( [ s_opt , sign_det/prod(s_opt) ] ) * V.';

  
  function [e,d] = ener( ss )
    
    p = prod(ss);
    p = sign_det/p;
    d = ( [ ss , p ] - s );
    
    e = d*d.';
    
    if nargout > 1
      
      d = d(1:end-1) - d(end) *( p./ss );

    end
    
  end
  
  
end
