function d2 = d24o(X, h, dim)
% Fourth-order accurate finite difference for SECOND derivative
% Same size as X
%
% X   : input array
% h   : grid spacing
% dim : differentiation dimension (1 or 2)

if nargin < 3
    dim = 1;
end

if size(X,dim) < 5
    error('d24o requires at least 5 points along dim=%d (got %d)', ...
          dim, size(X,dim));
end

d2 = zeros(size(X));

switch dim
    case 1  % along rows
        % Interior
        d2(3:end-2,:) = ( ...
            -1/12*X(1:end-4,:) ...
            +4/3 *X(2:end-3,:) ...
            -5/2 *X(3:end-2,:) ...
            +4/3 *X(4:end-1,:) ...
            -1/12*X(5:end,:) ) / h^2;

        % Forward boundary
        d2(1,:) = ( 35/12*X(1,:) - 26/3*X(2,:) + 19/2*X(3,:) ...
                   -14/3*X(4,:) + 11/12*X(5,:) ) / h^2;

        d2(2,:) = ( 11/12*X(1,:) - 5/3*X(2,:) + 1/2*X(3,:) ...
                    + 1/3*X(4,:) - 1/12*X(5,:) ) / h^2;

        % Backward boundary
        d2(end-1,:) = ( 11/12*X(end,:) - 5/3*X(end-1,:) + 1/2*X(end-2,:) ...
                         + 1/3*X(end-3,:) - 1/12*X(end-4,:) ) / h^2;

        d2(end,:) = ( 35/12*X(end,:) - 26/3*X(end-1,:) + 19/2*X(end-2,:) ...
                      -14/3*X(end-3,:) + 11/12*X(end-4,:) ) / h^2;

    case 2  % along columns
        % Interior
        d2(:,3:end-2) = ( ...
            -1/12*X(:,1:end-4) ...
            +4/3 *X(:,2:end-3) ...
            -5/2 *X(:,3:end-2) ...
            +4/3 *X(:,4:end-1) ...
            -1/12*X(:,5:end) ) / h^2;

        % Forward boundary
        d2(:,1) = ( 35/12*X(:,1) - 26/3*X(:,2) + 19/2*X(:,3) ...
                   -14/3*X(:,4) + 11/12*X(:,5) ) / h^2;

        d2(:,2) = ( 11/12*X(:,1) - 5/3*X(:,2) + 1/2*X(:,3) ...
                    + 1/3*X(:,4) - 1/12*X(:,5) ) / h^2;

        % Backward boundary
        d2(:,end-1) = ( 11/12*X(:,end) - 5/3*X(:,end-1) + 1/2*X(:,end-2) ...
                         + 1/3*X(:,end-3) - 1/12*X(:,end-4) ) / h^2;

        d2(:,end) = ( 35/12*X(:,end) - 26/3*X(:,end-1) + 19/2*X(:,end-2) ...
                      -14/3*X(:,end-3) + 11/12*X(:,end-4) ) / h^2;

    otherwise
        error('dim must be 1 or 2');
end
end