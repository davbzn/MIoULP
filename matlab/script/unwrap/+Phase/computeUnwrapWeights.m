function [w_Deltax w_Deltay] = computeUnwrapWeights(w)
ws = w.^2;
w_Deltax = min(ws(:, 1:end-1), ws(:, 2:end));
w_Deltay = min(ws(1:end-1, :), ws(2:end, :));
