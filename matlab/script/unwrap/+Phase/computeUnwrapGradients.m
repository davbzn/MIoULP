function [Deltax Deltay] = computeUnwrapGradients(psi)
Deltax = wrapToPi(psi(:, 2:end) - psi(:, 1:end-1));
Deltay = wrapToPi(psi(2:end, :) - psi(1:end-1, :));
