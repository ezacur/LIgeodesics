function R = projectOnSO( R )

  [u,s,v] = svd( R );
  R = u*diag( [ ones(1,size(R,1)-1) , det( u*v.' ) ] )*v.';

end
