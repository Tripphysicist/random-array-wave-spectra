function DS = imlm(Sxy,Gwt,thetai,fi,k,opt)
%IMLM  Iterated Maximum Likelihood Method based on DIWASP (Johnson 2002)
%
% CALL: DS = imlm(Sxy, Gwt, thetai, fi, k, opt);
%
%  DS     = Directional distribution (spreading function) size nt x nf
%  Sxy    = matrix of cross spectral densities size m x m x nf
%  Gwt    = matrix of transfer function (abs(Gwt)==1) size m x nt x nf
%  thetai = angle vector length nt
%  fi     = frequency vector length nf
%  k      = index vector to frequencies where Sf>0 length <= nf
%  opt    = options structure (maxiter, relax, errortol, etc.)
%
% -------------------------------------------------------------------------
% NOTE / PROVENANCE:
% This implementation replaces the legacy WAFO power-law formulation with the
% linear damped error-correction algorithm from DIWASP (Directional Wave
% Analysis Software Package, Johnson 2002; Krogstad 1988), as utilized in:
%
%   Simanesew, A. W., H. E. Krogstad, K. Trulsen, and J. C. Nieto Borge, 2018:
%   "Bimodality of directional distributions in ocean wave spectra: a comparison
%   of data-adaptive estimation techniques." J. Atmos. Oceanic Technol., 35 (2), 365-384.
%
% Algorithmic formulation:
%   1. Initial estimate: DS0 = MLM(Sxy)
%   2. For each iteration iz = 1:miter:
%      a. Re-synthesize cross-spectra Sxy_synth from current directional estimate DS
%      b. Compute MLM estimate T from Sxy_synth: T = MLM(Sxy_synth)
%      c. Apply linear update with momentum damping:
%           ei = gamma * ((DS0 - T) + alpha * (T - Told))
%           DS = DS + ei
%      d. Enforce non-negativity and normalize: int DS(theta, f) dtheta = 1
%
% Default parameters (matching DIWASP / Simanesew et al. 2018):
%   gamma    = 0.1   (relaxation factor / learning rate)
%   alpha    = 0.1   (momentum / inertial damping factor to suppress ringing)
%   miter    = 10    (number of iterations, preventing noise over-fitting)
%   errorTol = 1e-3  (convergence tolerance)
%
% The original WAFO implementation (undamped, quadratic power-law with Li=1.4,
% Bi=2.0) has been preserved under 'imlmw.m' for benchmarking.
% -------------------------------------------------------------------------

[m, nt, nf] = size(Gwt);

% Default algorithm parameters from DIWASP / Simanesew et al. (2018)
gamma    = 0.1;
alpha    = 0.1;
miter    = 10; % DIWASP standard
errorTol = 1e-3;
display  = 0;

if nargin > 5 && isstruct(opt)
    if isfield(opt, 'relax') && ~isempty(opt.relax),       gamma    = opt.relax;   end
    if isfield(opt, 'imlm_iter') && ~isempty(opt.imlm_iter), miter  = opt.imlm_iter; end
    if isfield(opt, 'errortol') && ~isempty(opt.errortol), errorTol = opt.errortol; end
    if isfield(opt, 'message') && ~isempty(opt.message),   display  = (opt.message > 1); end
end

% 1. Initial estimate from standard MLM
DS0 = mlm(Sxy, Gwt, thetai, fi, k, opt);
DS  = DS0;
T   = DS0;
Told = T;
DSold = DS;

Sxy_synth = zeros(m, m, nf);

% 2. Run DIWASP Iteration across active frequency bins
for iz = 1:miter
    % a. Re-synthesize cross-spectra from current estimate DS
    for ix = 1:m
        Sxy_synth(ix, ix, :) = trapz(thetai, squeeze(Gwt(ix,:,:).*conj(Gwt(ix,:,:))).*DS);
        for iy = (ix+1):m
            Sxy_synth(ix, iy, :) = trapz(thetai, squeeze(Gwt(ix,:,:).*conj(Gwt(iy,:,:))).*DS);
            Sxy_synth(iy, ix, :) = conj(Sxy_synth(ix, iy, :));
        end
    end
    
    % b. Re-estimate MLM from synthesized cross-spectra
    T = mlm(Sxy_synth, Gwt, thetai, fi, k, opt);
    
    % c. DIWASP update rule with momentum damping
    ei = gamma * ((DS0 - T) + alpha * (T - Told));
    Told = T;
    DS = DS + ei;
    
    % d. Enforce non-negativity and unit area normalization
    DS = normspfn(DS, thetai);
    
    % e. Check convergence
    tol = max(abs(DS(:) - DSold(:)));
    if display
        disp(['DIWASP IMLM Iteration ' num2str(iz) ' of ' num2str(miter) ' Error = ' num2str(tol)]);
    end
    if tol < errorTol
        break;
    end
    DSold = DS;
end

DS = normspfn(DS, thetai);
return; % imlm
