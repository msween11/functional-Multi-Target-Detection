function Estimate = get_phase_from_bispectrum_gap(B_mat,d)
% Given the bispectrum matrix and signal length, this function generates
% the phase of the original signal using spectral method with largest
% spectral gap described in Algorithm 1 of the paper
% The function asserts Hermitian by forcing the first row, first column
% and the diagonal to be 1.
%
% Jan 2018
% Hua Chen
% https://github.com/ARKEYTECT/Bispectrum_Inversion

        % assert Hermitian 
        B_mat = (B_mat+B_mat')/2;
        B_mat(:, 1) = 1;
        B_mat(1, :) = 1;
        B_mat = B_mat - diag(diag(B_mat)) + diag(ones(d, 1));
        [U,D] = eig(B_mat);
        eigval = real(diag(D));
        [eigval,sorted_id] = sort(eigval, 'descend');
        U = U(:, sorted_id);
        eigval = eigval.';
        dif = min([abs(diff(eigval)) max(eigval);max(eigval) abs(diff(eigval))]);
        [~,I] = sort(dif,'descend'); 
        V = U(:,I(1));
        V = V/V(1);
        Estimate = V./abs(V); 
        
end 