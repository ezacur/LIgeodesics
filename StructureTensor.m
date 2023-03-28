function C = StructureTensor( B )

  CHECK_BASIS = true;
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
  
  
  if 0 && CHECK_BASIS
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

  pB = pinv(B);

  tosq = @(x) reshape(x,N,N);
  C = zeros( K , K , K );
  for a=1:K
    for b=1:K
      
      COMM = vec( tosq( B(:,a) )*tosq( B(:,b) ) - tosq( B(:,b) )*tosq( B(:,a) ) );
      
      C(a,b,:) = pB*COMM(:);
     
    end
  end
  

if 0
  E = zeros(K,K);
  for a=1:K
    for b=1:K

E(a,b) = maxnorm( commutator( tosquare(B(:,a)),tosquare(B(:,b)) ) , tosquare( B*vec( C(a,b,:) ) ) );

    end
  end
end




end
