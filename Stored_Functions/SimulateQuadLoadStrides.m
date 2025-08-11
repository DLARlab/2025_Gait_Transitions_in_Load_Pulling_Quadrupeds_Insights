function [R, R9, T, Y, GRF, F_Load, P, X_accum_tru, num_strides] = SimulateQuadLoadStrides(X_accum)
%SIMULATEQUADLOADSTRIDES   simulate 1–N strides of the quad-load gait.
%
% [T,Y,GRF,F_Load,P] = simulateQuadLoadStrides(X_accum, num_strides)
%
%   X_accum      – length-44 + (num_strides–1)*4 vector
%   num_strides  – integer ≥1
%
% Outputs:
%   T       – stacked time
%   Y       – stacked states
%   GRF     – stacked ground reaction forces
%   F_Load  – stacked leash forces (empty if single-stride quad only)
%   P       – the accumulated P vector

  % fixed vs. extra per stride
  fixed_length = 44; varying_length = 13;


  num_strides = (length(X_accum) - fixed_length)/varying_length + 1;

  R = []; R9 = []; T = []; Y = []; GRF = []; F_Load = []; P = [];
  X_curr = X_accum(1:fixed_length);
  
  X_accum_tru = X_accum;

  for stride_idx = 1:num_strides
    % call the v2 function
    % X_curr(1:22) = EventTimingRegulation(X_curr(1:22));
    [R_curr, T_curr, Y_curr, P_curr, GRF_curr, F_curr] = Quad_Load_ZeroFun_Transition_v2(X_curr,'skipSolve');

    if stride_idx ==1
       X_accum_tru(14:22) = P_curr(1:9)';
    else
       X_accum_tru(fixed_length+1:fixed_length+9) =  P_curr(1:9)';
    end

    
    R = [R;R_curr];
    R9 = [R9;R_curr(1:9)];
    % stitch
    if ~isempty(T)
      T_curr = T_curr + T(end);
    end
    T    = [T;   T_curr];
    Y    = [Y;   Y_curr];
    GRF  = [GRF; GRF_curr];
    if ~isempty(F_curr)
      F_Load = [F_Load; F_curr];
      P      = [P; P_curr];  % last one only
    end

    if stride_idx < num_strides
        Xn = zeros(fixed_length, 1);
        Xn(1:13) = Y_curr(end, 2:14)';
        Xn(37) = Y_curr(end, 15) - Y_curr(end, 1);
        Xn(38) = Y_curr(end, 16);
        Xn(14:22) = X_accum(fixed_length + (stride_idx-1)*varying_length + 1 : fixed_length + (stride_idx-1)*varying_length + 9);
        Xn(24:27) = X_curr(28:31);
        Xn(28:31) = X_accum(fixed_length + (stride_idx-1)*varying_length + 10 : fixed_length + (stride_idx-1)*varying_length + 13);
        Xn([23 32:36 39:44]) = X_accum([23 32:36 39:44]);
        X_curr = Xn;
    end
  end
end
