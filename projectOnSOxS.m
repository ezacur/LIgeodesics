function R = projectOnSOxS( R )

  [u,s,v] = svd( R );
  R = u*diag( [ ones(1,size(R,1)-1) , det( u*v.' ) ] )*v.';
  
  R = R*trace(s)/size(s,1);

end
