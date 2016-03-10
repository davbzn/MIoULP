function c = weightedLaplacianFromGrad(Deltax, w_Deltax, Deltay,  w_Deltay)
% Matrix of weights of the Laplacian used in Ghiglia and Romero.
% Computes eq. (34) in Ghiglia and Romero JOSA A 11(1) 1994 p107, specified
% in terms of the phase differences Delta_x, Delta_y and their weights.
% Deltax is M x (N-1) and Deltay is (M-1) x N.
c = zeros(size(Deltax)+[0 1]);
c(:, 1:end-1) = w_Deltax .* Deltax;
c(:, 2:end) = c(:, 2:end) - w_Deltax .* Deltax;
c(1:end-1, :) = c(1:end-1, :) + w_Deltay .* Deltay;
c(2:end, :) = c(2:end, :) - w_Deltay .* Deltay;
