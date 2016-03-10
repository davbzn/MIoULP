function [phi r k] = reconFromGradLS2DPCD(Deltax,  w_Deltax, Deltay, w_Deltay, phi, opts)
% Reconstruct 2D phase from gradient, least squares, PCD method
%
%   PCD = preconditioned conjugate gradient
%   phi = reconFromGradLS2DPCD(Deltax,  w_Deltax, Deltay, w_Deltay) performs
%
%	Algorithm 3 of Ghiglia and Romero JOSA A 11(1) 1994 p107.
%   Deltax is M x (N-1) and Deltay is (M-1) x N. Deltax and Deltay are the
%   phase gradients. For unwrapping, we just use (2) and (3) for these. For
%   reconstruction, use the inverse of the variance of the measurement
%   error. The opts structure contains the options, with fields
%       max_iterations          [1000]
%       min_normr_ratio         [0.01]

% wmax = max(max(max(w_Deltax)), max(max(w_Deltay)));
% w_Deltax = w_Deltax / wmax;
% w_Deltay = w_Deltay / wmax;

if ~exist('phi', 'var') || isempty(phi)
    phi = zeros(size(Deltax) + [0 1]);
end
if ~exist('opts', 'var') 
    opts = struct;
end
opts = struct_defaults(opts, 'max_iterations', 1000, 'min_normr_ratio', 1e-6, 'monitor', false);
   

% Allow use of an initial guess
r = Phase.weightedLaplacianFromGrad(Deltax,  w_Deltax, Deltay, w_Deltay) - Phase.weightedLaplacianFromUnwrapped(phi, w_Deltax, w_Deltay);
r0 = sum(sum(abs_sqd(r)));

k = 0;
%while (k == 0 || max(max(abs(phim1 - phi)))) > opts.min_diff && k < opts.max_iterations 
normr_ratio = sum(sum(abs_sqd(r)))/r0;
while normr_ratio > opts.min_normr_ratio && k < opts.max_iterations 
    z = PoissonsViaDCT(r);
    
    % Step 3
    k = k + 1;
    if k > 1
        zm2 = zm1;
        rm2 = rm1;
        pm1 = p;
    end
    zm1 = z;
    rm1 = r;
    phim1 = phi;
    
    rm1_dot_zm1 = sum(sum(rm1 .* zm1));
    if k == 1
        % Step 4
        p = zm1;
    elseif k > 1
        % Step 5
        beta = rm1_dot_zm1/sum(sum(rm2 .* zm2));
        p = zm1 + beta * pm1;
    end
    % Step 6
    Qp = Phase.weightedLaplacianFromUnwrapped(p, w_Deltax, w_Deltay);
    alpha = rm1_dot_zm1 / sum(sum(p .* Qp));
    phi = phim1 + alpha*p;
    r = rm1 - alpha*Qp;
    normr_ratio = sum(sum(abs_sqd(r)))/r0;
    if opts.monitor
        sprintf('%.3f ', normr_ratio)
    end        
end
if opts.monitor
    sprintf('\n')
end
if k == opts.max_iterations
    warning('Stopped after maximum iterations');
end
