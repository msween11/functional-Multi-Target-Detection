function d1 = d14o(X, h, dim)
% Fourth-order accurate finite difference for first derivative
% Same size as X (like gradient)
%
% X   : input matrix
% h   : grid spacing
% dim : dimension (1 = rows, 2 = columns)

if nargin < 3
    dim = 1;
end

d1 = zeros(size(X));

switch dim
    case 1  % derivative along rows
        % Interior: centered 4th order
        d1(3:end-2,:) = ( ...
            (1/12)*X(1:end-4,:) ...
          - (2/3) *X(2:end-3,:) ...
          + (2/3) *X(4:end-1,:) ...
          - (1/12)*X(5:end,:) ) / h;

        % Left boundary (forward, 4th order)
        d1(1,:) = (-25/12*X(1,:) + 4*X(2,:) - 3*X(3,:) + 4/3*X(4,:) - 1/4*X(5,:)) / h;
        d1(2,:) = (-1/4*X(1,:) - 5/6*X(2,:) + 3/2*X(3,:) - 1/2*X(4,:) + 1/12*X(5,:)) / h;

        % Right boundary (backward, 4th order)
        d1(end,:)   = ( 25/12*X(end,:) - 4*X(end-1,:) + 3*X(end-2,:) ...
                        - 4/3*X(end-3,:) + 1/4*X(end-4,:) ) / h;
        d1(end-1,:) = ( 1/4*X(end,:) + 5/6*X(end-1,:) - 3/2*X(end-2,:) ...
                        + 1/2*X(end-3,:) - 1/12*X(end-4,:) ) / h;

    case 2  % derivative along columns
        % Interior
        d1(:,3:end-2) = ( ...
            (1/12)*X(:,1:end-4) ...
          - (2/3) *X(:,2:end-3) ...
          + (2/3) *X(:,4:end-1) ...
          - (1/12)*X(:,5:end) ) / h;

        % Bottom boundary (forward)
        d1(:,1) = (-25/12*X(:,1) + 4*X(:,2) - 3*X(:,3) + 4/3*X(:,4) - 1/4*X(:,5)) / h;
        d1(:,2) = (-1/4*X(:,1) - 5/6*X(:,2) + 3/2*X(:,3) - 1/2*X(:,4) + 1/12*X(:,5)) / h;

        % Top boundary (backward)
        d1(:,end)   = ( 25/12*X(:,end) - 4*X(:,end-1) + 3*X(:,end-2) ...
                        - 4/3*X(:,end-3) + 1/4*X(:,end-4) ) / h;
        d1(:,end-1) = ( 1/4*X(:,end) + 5/6*X(:,end-1) - 3/2*X(:,end-2) ...
                        + 1/2*X(:,end-3) - 1/12*X(:,end-4) ) / h;

    otherwise
        error('Only dim = 1 or 2 supported.');
end
end