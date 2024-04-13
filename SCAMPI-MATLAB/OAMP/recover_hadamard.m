function W_r = recover_hadamard (rp, N, Q)
H = hadamard(N);
for i = 1:1:Q
  W_r(i,:) = sqrt(1/N)*H(rp{1,1}(i),:);
end
