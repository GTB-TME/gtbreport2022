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
% different categories, e.g. infectious TB (see 'Read me' for more 
% detail on compartment lookups)
% - gps: lists nested stratifications, e.g. public/private sectors

% Outputs:
% ========
% - M: matrix.based model specification, with the following components:
% - M.lin captures linear transitions in the model
% - M.Dxlin captures linear transitions related to diagnosis and treatment
% initiation (used for scaling these rates in a time-dependent way to
% reflect COVID disruptions)
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
% - get_objective.m


function M = make_model(p, r, i, s, gps)

% --- Get the linear rates ------------------------------------------------
m    = zeros(i.nstates);
m2pu = zeros(i.nstates);
m2pr = zeros(i.nstates);

Lf   = i.Lf;
Ls   = i.Ls;
I    = i.I;
Txs  = i.Tx;
Rlo  = i.Rlo;
Rhi  = i.Rhi;
R    = i.R;

% --- Fast progression and LTBI stabilisation
source  = Lf; 
destins =             [I,            Ls]; 
rates   = [r.progression, r.LTBI_stabil];
m(destins, source) = m(destins, source) + rates';

% --- Reactivation
source = Ls; destin = I; rate = r.reactivation;
m(destin, source) = m(destin, source) + rate;

% --- Primary careseeking
source = I;
destin = Txs.pu;
rate   = r.Tx_init(1) + r.g;
m2pu(destin, source) = m2pu(destin, source) + rate;

source = I;
destin = Txs.pr;
rate   = r.Tx_init(2);
m2pr(destin, source) = m2pr(destin, source) + rate;


for ip = 1:length(gps.provs)
    prov = gps.provs{ip};
    Tx = Txs.(prov);
    
    % --- Treatment outcomes
    source  = Tx;
    destins = [Rlo               R                     Rhi];
    rates   = [r.Tx*p.cure(ip),  r.Tx*(1-p.cure(ip)),  r.default(ip)];
    m(destins, source) = m(destins, source) + rates';    
end

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
sources = s.infectious;
destin  = Rhi;
rates   = r.self_cure;
m(destin, sources) = m(destin, sources) + rates;

M.lin      = sparse(m - diag(sum(m,1)));
M.Dxlin_pu = sparse(m2pu - diag(sum(m2pu,1)));
M.Dxlin_pr = sparse(m2pr - diag(sum(m2pr,1)));


% --- Get the nonlinear rates ---------------------------------------------

% --- Allocating transitions
m = zeros(i.nstates);

m(i.Lf, [s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R]) = 1;
m(:,[s.Lf, s.Ls, s.Rlo, s.Rhi, s.R]) = m(:,[s.Lf, s.Ls, s.Rlo, s.Rhi, s.R])*(1-p.imm);
M.nlin = sparse(m - diag(sum(m,1)));


% --- Getting force-of-infection
m = zeros(1,i.nstates);
m(s.infectious) = r.beta;
m(:,setdiff(s.infectious,s.I)) = m(:,setdiff(s.infectious,s.I));
M.lambda = sparse(m);


% --- Get the mortality rates
m = zeros(i.nstates,2);
m(:,1) = r.mort;
m(s.infectious,2) = r.mort_TB;
% m([s.infectious, intersect(s.Tx, s.pr)]) = r.mort_TB;
M.mortvec = m;