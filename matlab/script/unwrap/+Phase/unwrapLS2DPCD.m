function [phi r k] = unwrapLS2DPCD(psi, w, phi, opts)
% Unwrap 2D phase using weighted preconditioned conjugate gradient method.
%
%	phi = unwrapLS2DPCD(psi, w, phi, opts)
%
%	performs Algorithm 3 of Ghiglia
%   and Romero JOSA A 11(1) 1994 p107. psi is the wrapped phase, and w is
%   the corresponding weights. phi is an initially guess, and can be omitted
%   left empty. opts is an options structure, defined in the
%   documentation for Phase.reconFromGradLS2DPCD.
%   the corresponding weights. It has fields
%       max_iterations          [1000]
%       min_normr_ratio         [0.01]
%
%	See also Phase.unwrapLS2DFourier Phase.unwrapLS2DPicard
%	Phase.reconFromGradLS2DPCD
if ~exist('w', 'var') w = ones(size(psi)); end
w = w / max(max(w));
if ~exist('phi', 'var') || isempty(phi)
    phi = zeros(size(psi));
end
if ~exist('opts', 'var') 
    opts = struct;
end

[Deltax Deltay] = Phase.computeUnwrapGradients(psi);
[w_Deltax w_Deltay] = Phase.computeUnwrapWeights(w);

[phi r k] = Phase.reconFromGradLS2DPCD(Deltax, w_Deltax, Deltay, w_Deltay, phi, opts);

% Add arbitrary constant for best agreement 
phi = phi - angle(mean(mean(w .* exp(1i*(phi - psi)))));
    
    
        



