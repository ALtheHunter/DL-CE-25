function generate_matrix(M, N)
G = rand(N,N) > 0.1;
H = hadamard(N);
% G = zeros(M,N);
% G = ones(N,N);
% G = H - G.* H;
G = G.* H;
save G;
%save(G, G, '-v7.3');
% X(randP) = X;
% W = G * X;
% save W;
disp('!');
end