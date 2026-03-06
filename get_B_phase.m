    function B_phase = get_B_phase(B_mat)
        
%         B_phase = B_mat./abs(B_mat);
%         B_phase(abs(B_mat)<1*10^-12) = 1;
          B_phase = exp(sqrt(-1)*angle(B_mat));
        
end