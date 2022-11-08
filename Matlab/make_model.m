% Given all model parameters, construct the full set of model equations in
% matrix form. This matrix-based model specification is then fed to the
% differential equation solver as an input.

% Arguments:
% ==========
% - p: Model parameters that are proportions
% - r: Model parameters that are per-capita rates
% - i: Compartment lookup giving index of each model compartment (see 'Read me'
% for more detail on compartment lookups)
% - s: Compartment lookup giving sets of compartment indices relating to
% different categories, e.g. HIV+ve vs HIV-ve (see 'Read me' for more 
% detail on compartment lookups)
% - gps: lists nested stratifications, e.g. HIV status, and public/private sectors
% - opts: opts.hiv set to 1 for countries with >10% of TB being HIV
% coinfected (0 otherwise) and opts.provs set to 1 for countries with a
% strong role for the private sector (0 otherwise)

% Outputs:
% ========
% - M: matrix.based model specification, with the following components:
% - M.lin captures linear transitions in the model
% - M.Dxlin captures linear transitions related to diagnosis and treatment
% initiation (used for scaling these rates in a time-dependent way to
% reflect COVID disruptions)
% - M.linHIV captures rates corresponding to HIV acquisition (used for
% capturing the time-dependent dynamics of the HIV epidemic)
% - M.nlin specifies source and origin compartments for all transitions
% arising from infection (e.g. from uninfected to 'latent fast'
% compartments)
% - M.lambda calculates the force-of-infection (equivalently, the ARTI)
% - M.mortvec captures all mortality rates 

% Dependencies:
% =============
% None

% This function is used by:
% =========================
% - BRR_Get_outputs.m
% - get_objective2D.m


function M = make_model(p, r, i, s, gps, opts)

% --- Get the linear rates ------------------------------------------------
m  = zeros(i.nstates);                                                     % Matrix for linear transitions in the model (i.e. except infection)
m2 = zeros(i.nstates);                                                     % Separate linear matrix for linkage to treatment - this will be multiplied by a time-dependent factor in 'goveqs_basis', to reflect COVID-related disruptions

for ih = 1:length(gps.hiv)
    
    hiv = gps.hiv{ih};
    getst = @(st) i.(st).(hiv);
    
    Lf   = getst('Lf');
    Ls   = getst('Ls');
    I    = getst('I');
    Txs  = getst('Tx');
    Rlo  = getst('Rlo');
    Rhi  = getst('Rhi');
    R    = getst('R');
    
    % --- Fast progression and LTBI stabilisation
    source  = Lf;
    destins = [I,                 Ls];
    rates   = [r.progression(ih), r.LTBI_stabil(ih)];
    m(destins, source) = m(destins, source) + rates';
    
    % --- Reactivation
    source = Ls; destin = I; rate = r.reactivation(ih);
    m(destin, source) = m(destin, source) + rate;
    
    % --- Primary careseeking
    source = I;
    destin = Txs.pu;
    rate   = r.Tx_init(1);
    m2(destin, source) = m2(destin, source) + rate;
    
    source = I;
    destin = Txs.pr;
    rate   = r.Tx_init(2);
    m2(destin, source) = m2(destin, source) + rate;

    % --- Suppleemntary case-finding (assumed from 2014 onwards)
    source = I;
    destin = Txs.pu;
    rate   = r.casefinding;
    m2(destin, source) = m2(destin, source) + rate;
    
    % --- Treatment outcomes
    source  = Txs.pu;
    destins = [Rlo    Rhi];
    rates   = [r.Tx,  r.default];
    m(destins, source) = m(destins, source) + rates';

    source  = Txs.pr;
    destins = [Rlo    Rhi];
    rates   = [r.Tx,  r.default];
    m(destins, source) = m(destins, source) + rates';

    % --- Relapse
    sources = [Rlo, Rhi, R];
    destin  = I;
    rates   = r.relapse;
    m(destin, sources) = m(destin, sources) + rates;
    
    sources = [Rlo, Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;
    
    % --- Self cure
    sources = intersect(s.infectious,s.(hiv));
    destin  = Rhi;
    rates   = r.self_cure(ih);
    m(destin, sources) = m(destin, sources) + rates;
    
end

m3 = zeros(i.nstates);                                                     % Separate linear matrix for HIV incidence - this will be multiplied by a time-dependent factor in 'goveqs_basis', to reflect HIV acquisition
if opts.hiv
    % --- HIV acquisition                                                  
    sources = s.h0;
    destins = s.h1;
    inds    = sub2ind(size(m), destins, sources);
    rates   = 1;
    m3(inds) = m3(inds) + rates;
    
    % --- ART initiation
    sources = s.h1;
    destins = s.hart;
    inds    = sub2ind(size(m), destins, sources);
    rates   = r.ART_init;
    m(inds) = m(inds) + rates;
end

% --- Bring them together
M.lin    = sparse(m - diag(sum(m,1)));
M.Dxlin  = sparse(m2 - diag(sum(m2,1)));
M.linHIV = sparse(m3 - diag(sum(m3,1)));                       


% --- Get the nonlinear rates ---------------------------------------------

% --- Allocating transitions
m = zeros(i.nstates);
for ih = 1:length(gps.hiv)
    hiv = gps.hiv{ih};
    sources = intersect([s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(hiv));
    destin  = i.Lf.(hiv);
    m(destin, sources) = m(destin, sources) + 1;
    % Adjust for any immune protection
    cols = intersect([s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(hiv));
    m(:,cols) = m(:,cols)*(1-p.imm(ih));
end
% Relative acquisition for HIV
cols = [s.h1, s.hart];
m(:,cols) = m(:,cols)*p.HIVlam;
M.nlin = sparse(m - diag(sum(m,1)));


% --- Getting force-of-infection
m = zeros(1,i.nstates);
m(intersect(s.infectious,s.h0)) = r.beta(1);
m(intersect(s.infectious,[s.h1,s.hart])) = r.beta(2);
M.lambda = sparse(m);


% --- Get the mortality rates
m         = zeros(i.nstates,3);

% First column is for non-TB-related mortality
m(:,1)    = r.mort;                                                         
m(s.h1,1) = r.HIV_mort;                                                    

% Second column is for HIV-ve TB mortality
inds = intersect(s.h0, s.infectious);
m(inds,2) = r.mort_TB(1);                                          

% Third column is for HIV+ve TB mortality
inds = intersect([s.h1,s.hart], s.infectious);
m(inds,3) = r.mort_TB(2);                                                  
M.mortvec = m;