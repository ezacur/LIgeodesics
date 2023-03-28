function [d,U] = Log_GLr( Q , varargin )

  for va=1:numel(varargin)
    if isnumeric(varargin{va}) && ~isscalar(varargin{va})
      varargin{va} = -varargin{va};
    end
  end

  Q = Q \ eye(size(Q));
  [d,U] = Log_GL( Q , varargin{:} );
  U = -U;

end
