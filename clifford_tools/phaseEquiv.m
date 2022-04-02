function output = phaseEquiv(U,V)
%%%% check two matrices are equivalent up to a phase
tolerance = 1e-7;
phase_tol = 1e-4;

[max_phase,index] = max(abs(U(:)));

[row,col] = ind2sub(size(U),index);

Uphase = U(row,col);

Vphase = V(row,col);

if or(abs(Uphase)<phase_tol,abs(Vphase)<phase_tol) ;
    output = false;
    return
end

U_corrected = (1/Uphase)*U;
V_corrected = (1/Vphase)*V;

diff = U_corrected - V_corrected;

sum_diffs = sum(sum(abs(diff)));

result = sum_diffs < tolerance;

output = result;

end
