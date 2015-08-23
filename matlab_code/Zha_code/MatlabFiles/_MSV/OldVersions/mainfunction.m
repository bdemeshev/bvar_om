% The constant paramater model written in the canonical form:
%  A x_t = B x_{t-1} + Gamma_u u_t + Pi pi_t
%   where f_t = [pi_t, y_t, R_t, mu_w_t, nu_t, E_t pi_{t+1}, E_t y_{t+1}]';
%         ut = [Sst*x_t{-1}; diag(e_t(1),...,e_t(h*))]
%         e_t is a 3-by-1 vector of fundamental shocks: e_{r t}, e_{wt}, and e_{nu t};
%         pi_t is a 2-by-1 vector of expectationsl errors for inflation and output: pi_{pi t} and pi_{y t}.
%
% Output format: x_t = G1*x_{t-1} + impact*e_t.
%
% Look for <<>> for ad hoc changes.
%See Notes pp.41-42.

smallval = 1.0e-8;  %<sqrt(machine epsilon).
locs_normalizedvars = [3 4 5]; %R_t, mu_w_t, nu_t;
seednumber = 7910;    %472534;   % if 0, random state at each clock time
if seednumber
   randn('state',seednumber);
   rand('state',seednumber);
else
   randn('state',fix(100*sum(clock)));
   rand('state',fix(100*sum(clock)));
end
n_normalizedvars = length(locs_normalizedvars);



%--- Reading in the values of the parameters.
%addpath D:\ZhaData\WorkDisk\LiuWZ\CalibrationMatlab\LWZmodel\Paravals2r
cd ./Paravals2r
%paravals2r_bench
%paravals2r_complex
paravals2r_noexistence
%paravals2r_try
cd ..

Q = [
 q11  1-q22
 1-q11  q22
     ];
disp(' ')
disp('Ergodic probability:')
qbar = fn_ergodp(Q)
g_etabar = qbar'*g_eta;
if (flag_swap_rs)
   g_eta = flipud(g_eta);
   g_gamma = flipud(g_gamma);
   g_rhor = flipud(g_rhor);
   g_phipi = flipud(g_phipi);
   g_phiy = flipud(g_phiy);
   Q = flipud(fliplr(Q));
end


%--- Composite regime-switching parameters.  The symbol "td" stands for "tilde."
htd = 1; %4;
Qtd = [
   Q(1,1) Q(1,1) 0       0
   0      0      Q(1,2)  Q(1,2)
   Q(2,1) Q(2,1) 0       0
   0      0      Q(2,2)  Q(2,2)
    ];
%--- psi1(std_t) = g_psi1(s_t,s_{t-1})
psi1td = [
   g_etabar*(1-g_eta(1)) / (g_eta(1)*(1-g_eta(1)))
   g_etabar*(1-g_eta(2)) / (g_eta(2)*(1-g_eta(1)))
   g_etabar*(1-g_eta(1)) / (g_eta(1)*(1-g_eta(2)))
   g_etabar*(1-g_eta(2)) / (g_eta(2)*(1-g_eta(2)))
   ];
%--- psi2(std_t) = g_psi2(s_{t-1})
psi2td = [
   (1-g_beta*g_etabar)*(1-g_eta(1)) / (g_eta(1)*(1+g_thetap*(1-g_alpha)/g_alpha));
   (1-g_beta*g_etabar)*(1-g_eta(2)) / (g_eta(2)*(1+g_thetap*(1-g_alpha)/g_alpha));
   (1-g_beta*g_etabar)*(1-g_eta(1)) / (g_eta(1)*(1+g_thetap*(1-g_alpha)/g_alpha));
   (1-g_beta*g_etabar)*(1-g_eta(2)) / (g_eta(2)*(1+g_thetap*(1-g_alpha)/g_alpha));
    ];
%--- gamma0(std_t) = g_gamma(s_{t})
gamma0td = [
   g_gamma(1)
   g_gamma(1)
   g_gamma(2)
   g_gamma(2)
    ];
%--- gamma1(std_t) = g_gamma(s_{t-1})
gamma1td = [
   g_gamma(1)
   g_gamma(2)
   g_gamma(1)
   g_gamma(2)
    ];
%--- rhor(std_t) = g_rhor(s_{t})
rhortd = [
   g_rhor(1)
   g_rhor(1)
   g_rhor(2)
   g_rhor(2)
    ];
%--- phipi(std_t) = g_phipi(s_{t})
phipitd = [
   g_phipi(1)
   g_phipi(1)
   g_phipi(2)
   g_phipi(2)
    ];
%--- phiy(std_t) = g_phiy(s_{t})
phiytd = [
   g_phiy(1)
   g_phiy(1)
   g_phiy(2)
   g_phiy(2)
    ];

%-------------------- Filling in A(std_t)'s, B(std_t)'s, etc for an expanded gensys. -------------------%
Astd = cell(htd,1);
Bstd = cell(htd,1);
Psistd = cell(htd,1);
for si=1:htd
   a1td = (1+g_beta*psi1td(si)*gamma0td(si));
   a2td = psi2td(si)*((g_xi+1)/g_alpha + b/(g_lambda-b));
   Astd{si} = [
     -a1td                        a2td                       0                      psi2td(si)  psi2td(si)*(b/(g_lambda-b))  g_beta*psi1td(si)      0
      0                          -(g_lambda+b)/g_lambda     -(g_lambda-b)/g_lambda  0           g_rhov - b/g_lambda          (g_lambda-b)/g_lambda  1
    -(1-rhortd(si))*phipitd(si)  -(1-rhortd(si))*phiytd(si)  1                      0           0                            0                      0
      0                           0                          0                      1           0                            0                      0
      0                           0                          0                      0           1                            0                      0
          ];
   Bstd{si} = [
     -gamma1td(si)  psi2td(si)*b/(g_lambda-b)  0           0       0       0     0
      0            -b/g_lambda                 0           0       0       0     0
      0             0                          rhortd(si)  0       0       0     0
      0             0                          0           g_rhow  0       0     0
      0             0                          0           0       g_rhov  0     0
          ];
   Psistd{si} = [
     0         0         0
     0         0         0
     g_sigmar  0         0
     0         g_sigmaw  0
     0         0         g_sigmav
       ];
end
%--- <<>>Reordering
if (flag_swap_fps)
   ki = 4;
   kk = 1;
   Atmp = Astd{ki};
   Astd{ki} = Astd{kk};
   Astd{kk} = Atmp;
   Btmp = Bstd{ki};
   Bstd{ki} = Bstd{kk};
   Bstd{kk} = Btmp;
   SwapM = zeros(ki,ki);
   SwapM(kk,ki) = kk;
   SwapM(ki,1) = 1;
   SwapM(2,2) = 1;
   SwapM(3,3) = 1;
   Qtd = SwapM*Qtd*SwapM;
end




Aind  = cell(htd,1);
Bind = cell(htd,1);
Psiind = cell(htd,1);
Piind = cell(htd,1);
G1ind = cell(htd,1);
impactind = cell(htd,1);
gevind = cell(htd,1);
gevind_abs = cell(htd,1);
euind = cell(htd,1);
for si=1:htd
   Aind{si} = [
      Astd{si}
      1 0 0 0 0 0 0
      0 1 0 0 0 0 0
      ];
   Bind{si} = [
      Bstd{si}
      0 0 0 0 0 1 0
      0 0 0 0 0 0 1
      ];
   Psiind{si} = [
      Psistd{si}
      0 0 0
      0 0 0
      ];
   Piind{si} = [
      zeros(5,2)
      1 0
      0 1
      ];

  disp('---------- Individual (possibly swapped) regimes in isolation --------------')
  disp('Under Regime')
  si
  
  %[G1ind{si},Cind,impactind{si},fmat,fwt,ywt,gevind{si},euind{si}]=gensys(Aind{si}, Bind{si}, zeros(7,1), Psiind{si}, Piind{si}, divind(si));
  %[G1ind{si},impactind{si},gevind{si},euind{si}]=msv_one(Aind{si}, Bind{si}, Psiind{si}, Piind{si});
  [G1ind{si},impactind{si},gevind{si},euind{si}]=msv_simple(Aind{si}, Bind{si}, Psiind{si}, Piind{si});
  [G1 impact gev eu]=msv_all(Aind{si}, Bind{si}, Psiind{si}, Piind{si});
  %[G1ind{si},impactind{si},v2junk,gevind{si},euind{si}]=gensys_dw(Aind{si}, Bind{si}, Psiind{si}, Piind{si}, divind(si));
  
  k=size(eu,2);
  for i=1:k
      gev_abs = abs(gev{i}(:,1).\gev{i}(:,2));
      gev_abs
      G1{i}
      impact{i}
      eu{i}
  end
end
% disp(' ')
% disp('----------------- Existence and uniqueness for each individual regime -----------------------')
% for si=1:htd
%    euind{si}
% end
return


