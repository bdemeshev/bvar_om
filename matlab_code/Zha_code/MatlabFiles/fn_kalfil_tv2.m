function [loglh,zt_tm1,Pt_tm1, loglh_t_allvalues] = fn_kalfil_tv2(Y_T, a, H, R, G, b, F, V, indxIni, indxDiffuse, z0,P0)
%[loglh,zt_tm1,Pt_tm1, loglh_t_allvalues] = fn_kalfil_tv(Y_T, a, H, R, G, b, F, V, indxIni, indxDiffuse, z0,P0)
%Time-varying Kalman filter (conditional on all the regimes in a Markov-switching model).  It computes
%  a sequence of one-step predictions and their covariance matrices, and the log likelihood.
%  The function uses a forward recursion algorithm.
%
% Note: for b, F, V, zt_tm1, and Pt_tm1, there is not need for the T+1 length because the values at T+1 have
%        never been used in the likelihood function.  See tz_kalfiltv() in kalman.c for an illustration.
%
%   State space model is defined as follows:
%       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)
%       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)
%     where a(t), H(t), b(t), and F(t) depend on s_t that follows a Markov-chain process and are taken as given.
%
%   Inputs are as follows:
%      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].
%        a is an n_y-by-T matrix of time-varying input vectors in the measurement equation.
%        H is an n_y-by-n_z-by-T 3-D of time-varying matrices in the measurement equation.
%        R is an n_y-by-n_y-by-T 3-D of time-varying covariance matrices for the error in the measurement equation.
%        G is an n_z-by-n_y-by-T 3-D of time-varying E(eta_t * eps_t').
%        ------
%        b is an n_z-by-(T+1) matrix of time-varying input vectors in the state equation with b(:,1) as an initial condition.
%        F is an n_z-by-n_z-by-(T+1) 3-D of time-varying transition matrices in the state equation with F(:,:,1) as an initial condition.
%        V is an n_z-by-n_z-by-(T+1) 3-D of time-varying covariance matrices for the error in the state equation with V(:,:,1) as an initial condition.
%        ------
%        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0
%                      (thus, indxDiffuse is automatically not enativated).;
%                 0: creating the initial condition depending on indxDiffuse (thus, indxDiffuse is active).
%        indxDiffuse: 1: using the diffuse initial condition by setting P=100*I and z=0;
%                     0: using the unconditional mean and variance for any given regime at time 0.
%        z0 is an n_z-by-1 vector of initial condition when indxIni=1. (Do not enter if indxIni=0.)
%        P0 is an n_z-by-n_z matrix of initial condition when indxIni=1. (Do not enter if indxIni=0.)
%
%   Outputs are as follows:
%      loglh is a value of the log likelihood function of the state-space model
%                                under the assumption that errors are multivariate Gaussian.
%      zt_tm1 is an n_z-by-(T+1) matrices of one-step predicted state vectors with z0_0m1 as a initial condition.
%      Pt_tm1 is an n_z-by-n_z-by-(T+1) 3-D of covariance matrices of zt_tm1 with P0_0m1 as a initial condition.
%      loglh_t_allvalues is a T-by-1 vector of loglh at time t for t=1:T.
%
%   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:
%             z0_0m1 = (I-F(:,:,1))\b(:,1)
%        vec(P0_0m1) = (I-kron(F(:,:,1),F(:,:,1)))\vec(V(:,:,1))
%   Note that all eigenvalues of the matrix F(:,:,1) are inside the unit circle when the state-space model is bounded (stationary).
%
%   March 2007, written by Tao Zha
%   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.
%==========================================================================
% Revision history:
%
%
%
%==========================================================================

% Copyright (C) 1997-2012 Tao Zha
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

Tp1 = size(Y_T,2) + 1;  %T+1;
n_y = size(a,1);
n_z = size(b,1);
%--- Checking input matrix dimensions
if (size(Y_T,1)~=n_y)
  error('kalf_tv(): Y_T and a must have the same number of rows')
end
%--- Allocating memory.
zt_tm1 = zeros(n_z,Tp1);
Pt_tm1 = zeros(n_z,n_z,Tp1);
loglh_t_allvalues = zeros(Tp1-1,1);
%--- Initializing.
loglh = 0.0;
if (indxIni)
   zt_tm1(:,1) = z0;
   Pt_tm1(:,:,1) = P0;
else
   if (indxDiffuse)
      %See Koopman and Durbin, "Filtering and Smoothing of State Vector for Diffuse State-Space Models," J. of Time Series Analysis, Vol 24(1), pp.85-99.
      zt_tm1(:,1) = zeros(n_z, 1);
      Pt_tm1(:,:,1) = 100*eye(n_z);
   else
      eigmax4F = max(abs(eig(F(:,:,1))));
      format long e
      eigmax4F
      if (eigmax4F < 1.0)
         zt_tm1(:,1) = (eye(n_z)-F(:,:,1))\b(:,1);
         V1 = V(:,:,1);
         Pt_tm1(:,:,1) = reshape((eye(n_z^2)-kron(F(:,:,1),F(:,:,1)))\V1(:) ,n_z,n_z);
      else
         %--- Do NOT use the following option.  It turns out that this will often generate explosive conditional liklihood
         %---   at the end of the sample, because Pt_tm1 shrinks to zero overtime due to the sigularity of the initila condition P_{1|0}.
         % zt_tm1(:,1) = zeros(size(zt_tm1(:,1)));
         % Pt_tm1(:,:,1) = V(:,:,1); %+eye(size(V(:,:,1)));

         nearinfinity = -1.0e+300
         loglh = nearinfinity;
         return; %Eearly exit.
         %error('kalf_tv(): For non-stationary solutions, the initial conditions must be supplied by, say, input arguments')
      end
   end
end


%====== See p.002 in LiuWZ. ======
indx_badlh = 0;   %1: bad likelihood with, say, -infinity of the LH value.
for t=2:Tp1
   tdata = t-1;

   %--- Setup.
   Htdata = H(:,:,tdata);
   Htdatatran = Htdata';
   ztdata = zt_tm1(:,tdata);
   Ptdata = Pt_tm1(:,:,tdata);
   PHtran_tdata = Ptdata*Htdatatran;
   Ft = F(:,:,t);
   Fttran = Ft';

   %--- Data.
   etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata;
   etdatatran = etdata';
   Dtdata = Htdata*PHtran_tdata + R(:,:,tdata);
   Dtdata = 0.5*(Dtdata + Dtdata');  %Making it symmetric against some rounding errors.
                      %This making-symmetric is very IMPORTANT; otherwise, we will get the matrix being singular message
                      %    and eigenvalues being negative for the SPD matrix, etc.  Then the likelihood becomes either
                      %    a bad number or a complex number.

   %--- State (updating).
   Kt_tdata = (Ft*PHtran_tdata+G(:,:,tdata))/Dtdata;
   Kt_tdatatran = Kt_tdata';
   zt_tm1(:,t) = b(:,t) + Ft*zt_tm1(:,tdata) + Kt_tdata*etdata;
   Pt_tm1(:,:,t) = Ft*Ptdata*Fttran - Kt_tdata*Dtdata*Kt_tdatatran + V(:,:,t);

   %--- Forming the log likelihood.
   detDtdata = det(Dtdata);
   %if (~isfinite(detDtdata))
   if (detDtdata < realmin)
      indx_badlh = 1;
      break;
   else
      loglh_tdata = -(0.5*n_y)*log(2*pi) - 0.5*log(detDtdata) - 0.5*(etdatatran/Dtdata)*etdata;
      loglh = loglh + loglh_tdata;
      loglh_t_allvalues(tdata) = loglh_tdata;
   end
end

if (indx_badlh)
   nearinfinity = -1.0e+300;
   loglh = nearinfinity;
   loglh_t_allvalues(tdata) = nearinfinity;
end

