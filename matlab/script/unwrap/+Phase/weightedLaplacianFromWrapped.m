function c = weightedLaplacianFromWrapped(psi, w_Deltax, w_Deltay)
% Computes eq. (34) in Ghiglia and Romero JOSA A 11(1) 1994 p107.
%   This uses the wrapped phase, and is intended for evaluation of the
%   constant term in the equations i.e. c.
%   See also Phase.weightedLaplacianFromGrad.
Deltax = wrapToPi(psi(:, 2:end) - psi(:, 1:end-1));
Deltay = wrapToPi(psi(2:end, :) - psi(1:end-1, :));
% ws = weights.^2;
% wxs = min(ws(:, 1:end-1), ws(:, 2:end));
% wys = min(ws(1:end-1, :), ws(2:end, :));
c = Phase.weightedLaplacianFromGrad(Deltax, w_Deltax, Deltay, w_Deltay);
