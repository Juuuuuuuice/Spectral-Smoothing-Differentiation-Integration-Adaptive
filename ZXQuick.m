function S = ZX(RAW, DA, D1, D2)

% ZX: Zone eXtraction method for adaptive spectral smoothing and calculus.
% A hierarchical framework for manifold decomposition and signal reconstruction.
%
% Copyright (c) 2026. All rights reserved. 
% Reference: https://doi.org/10.1016/j.measurement.2026.122204
%
% Description: 
% This script provides the implementation of the ZX method as described
% in the accompanying publication.

% Input:
%   RAW - Raw spectral data (column vector)
%   DA  - Order of calculus (0: smoothing, >0: derivative, <0: integral)
%   D1  - Global suppression parameter (Default: 2)
%   D2  - Local refinement parameter (Default: 1)


    %% --- 1. Initialization and Default Parameters ---
    RAW = RAW(:); 
    if nargin < 2, DA = 0; D1 = 2; D2 = 1; end 
    if nargin < 3, D1 = 2; D2 = 1; end 
    if nargin < 4, D2 = 1; end 

    %% --- 2. Hierarchical Processing ---
    % Stage A: Global Manifold Extraction
    frf_global = FindZone(RAW);
    frf_g_D1   = abs(frf_global .^ D1);        % [OPT] compute once, reuse below
    Sa_global  = TimeConv(RAW', frf_g_D1, 0)';

    % Stage B: Feature-Preserving Residual Hard-Thresholding
    DX  = RAW - Sa_global; 
    DDX = movstd(DX, [0, 10]); 
    is_noise   = DDX < (max(DDX) / 3);
    DX_refined = DX;
    DX_refined(is_noise) = 0;                   % keep DX intact for Sa2

    % Stage C: Secondary Sensing and Final Reconstruction
    frf_local = FindZone(DX_refined(3:end-2)); 
    Sa1 = TimeConv(RAW', frf_g_D1,                  DA);  % [OPT] reuse frf_g_D1
    Sa2 = TimeConv(DX',  abs(frf_local .^ D2),       DA);  
    
    S = Sa1' + Sa2';
end

%% ================= PROFESSIONAL SUB-MODULES =================

function cf = FindZone(Signal)
    [X, Y] = PickPoints(Signal);
    XY        = sortrows([X(:), Y(:)]);
    fit_model = CureFit(XY(1:end-6, 1), XY(1:end-6, 2));
    full_curve = fit_model(1:10000);
    cf = (full_curve(:) - full_curve(end)) ./ full_curve(:);
end

function [X, Y] = PickPoints(Signal)
    M = fft(Signal, 20001);
    Amplitude = abs(M); 
    A = Amplitude(1:10000) + flipud(Amplitude(10002:end));
    [~, mb] = min(movstd(A, 20));
    
    A_subset = A(1:mb);
    is_envelope_point = (A_subset >= flipud(cummax(flipud(A_subset))));
    
    X = find(is_envelope_point)'; 
    Y = A_subset(is_envelope_point)' / std(A_subset);
end

function S = TimeConv(Mraw, cf, DN)
    L  = 100;
    nc = size(Mraw, 2);
    Mr = [repmat(Mraw(:,1), 1, L), Mraw, repmat(Mraw(:,end), 1, L)]; % [OPT] repmat

    SC0 = GenKernel(cf(:)); 
    SC  = GLfd(SC0, DN);    
    SC1 = SC(5001 + floor(DN) : end - 5000);

    S0 = conv2(Mr, SC1(:)', 'same');
    S2 = S0(:, L+1 : L+nc);
    S1 = S0(:, L+2 : L+nc+1);
    ph = DN - floor(DN/2)*2;
    S  = ph*S1/2 + (2 - ph)*S2/2;
end

function TC = GenKernel(cf)
    IGG     = (cf(:) - min(cf)) / (max(cf) - min(cf));
    YT      = ifft(IGG, 20000);
    TC_vals = real(YT(1:10000));
    TC_vals = TC_vals - mean(TC_vals(9980:9999));
    TC = [flipud(TC_vals(2:9001)); TC_vals(1:9000)];
    TC = TC / sum(TC);
end

function dy = GLfd(y, v)
    n = length(y);
    k = (1:n-1)';
    w = [1; cumprod(1 - (v+1) ./ k)];           % [OPT] cumprod replaces for-loop
    res = conv(y(:), w, 'full');
    dy  = res(1:n)';
end

function fitr = CureFit(X, Y)
    % Fits the envelope curve using IRLS-LAR (no Curve Fitting Toolbox)
    % Model: a*exp(-(b*x)^p)+c
    % Parameters in alphabetical order: a(1), b(2), c(3), p(4)
    model = @(p, x) p(1) .* exp(-(p(2) .* x) .^ p(4)) + p(3);
    p0 = [4,   0.5, 0,     0.3 ];
    lb = [1,   0,   0.001, 0.25];
    ub = [20,  1,   2,     4   ];
    Xc = X(:);  Yc = Y(:);
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
                        'FunctionTolerance', 1e-6, 'MaxIterations', 400);

    % Warm start with plain OLS, then IRLS-LAR (mimics Robust='LAR')
    params = lsqcurvefit(model, p0, Xc, Yc, lb, ub, opts);
    for iter = 1:12
        resid    = Yc - model(params, Xc);
        w        = 1 ./ max(abs(resid), 0.5 * std(resid));
        sw       = sqrt(w / mean(w));
        p_new    = lsqcurvefit(@(p,x) model(p,x) .* sw, params, ...
                               Xc, Yc .* sw, lb, ub, opts);
        if norm(p_new - params) < 1e-6 * norm(params), break; end  % [OPT] early exit
        params   = p_new;
    end
    fitr = @(xq) model(params, xq(:));
end