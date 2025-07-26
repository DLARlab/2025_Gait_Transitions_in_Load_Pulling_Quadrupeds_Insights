function [total_cost, cost_terms, R2] = fms_NStridesOptimization_Quad_Load_v2(X_accum, t_exp, ft_exp, loading_force_exp, term_weights)

  num_strides = size(t_exp,1);
  fixed_length   = 44;    % still 44 per stride, since you're still carrying E(9) in X_multi
  varying_length  = 13;
  expected_len= fixed_length + max(0,num_strides-1)*varying_length;
  if numel(X_accum)~=expected_len
    error('bad X_multi length');
  end

  % split out the combined-vector arguments for the first stride
  X_curr = X_accum(1:fixed_length);

  % accumulators
  sd_diff_sum        = 0;
  ft_diff_sum        = 0;
  load_force_diff    = 0;
  

  SD_exp = [];
  SD_sim = [];
  FT_exp = [];
  FT_sim = [];
  LF_exp = [];
  LF_sim = [];

  for stride_idx = 1:num_strides
    
    [~, T_curr, Y_curr, P_curr, ~, F_leash_curr, ~] = Quad_Load_ZeroFun_Transition_v2(X_curr);

   
    % stride‐duration
    if iscell(t_exp)
        T_curr_dsample = DataLengthResampler(T_curr, t_exp{stride_idx,:});
        sd_diff_sum = sd_diff_sum + norm(t_exp{stride_idx,:} - T_curr_dsample);
        SD_exp = [SD_exp t_exp{stride_idx,:}];
        SD_sim = [SD_sim T_curr_dsample];
    else
        T_curr_dsample = DataLengthResampler(T_curr, t_exp(stride_idx,:));
        sd_diff_sum = sd_diff_sum + norm(t_exp(stride_idx,:) - T_curr_dsample);
        SD_exp = [SD_exp t_exp(stride_idx,:)];
        SD_sim = [SD_sim T_curr_dsample];
    end
    


    % foot‐touch times normalized by P(9)
    perm = [1 2 3 4 7 8 5 6];       % new order
    ft_sim = P_curr(perm) ./ P_curr(9);

    ft_dog = ft_exp(stride_idx, :);
    ft_diff_sum = ft_diff_sum + norm(ft_sim - ft_dog);
    FT_exp = [FT_exp ft_dog];
    FT_sim = [FT_sim ft_sim+(stride_idx-1)];


    % loading‐force
    if iscell(loading_force_exp)
      LFt = DataLengthResampler(loading_force_exp{stride_idx}, T_curr)';
    else
      LFt = DataLengthResampler(loading_force_exp(:,stride_idx), T_curr)';
    end
    load_force_diff = load_force_diff + norm(LFt - F_leash_curr);
    LF_exp = [LF_exp; LFt];
    LF_sim = [LF_sim; F_leash_curr];

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

  % package cost_terms
  cost_terms.strideduration = sd_diff_sum;
  cost_terms.ft             = ft_diff_sum;
  cost_terms.loadingforce   = load_force_diff;


  total_cost =...
    term_weights.strideduration * cost_terms.strideduration + ...
    term_weights.ft             * cost_terms.ft         + ...
    term_weights.loadingforce   * cost_terms.loadingforce;
    

    % Calculate the residual sum of squares (RSS) for Footfall Timings
    RSS_sd = sum((SD_exp - SD_sim).^2);
    TSS_sd = sum((SD_exp - mean(SD_exp)).^2);

    % Calculate the residual sum of squares (RSS) for Footfall Timings
    RSS_ft = sum((FT_exp - FT_sim).^2);
    TSS_ft = sum((FT_exp - mean(FT_exp)).^2);

    % Calculate the residual sum of squares (RSS) for Loading Force
    RSS_lf = sum((LF_exp - LF_sim).^2);
    TSS_lf = sum((LF_exp - mean(LF_exp)).^2);
    
    R2.strideduration = 1 - (RSS_sd / TSS_sd);
    R2.ft = 1 - (RSS_ft / TSS_ft);
    R2.loadingforce = 1 - (RSS_lf / TSS_lf);

    R2.weighted = (term_weights.strideduration * R2.strideduration + ...
                   term_weights.ft * R2.ft + ...
                   term_weights.loadingforce * R2.loadingforce)/(term_weights.strideduration + term_weights.ft + term_weights.loadingforce);
end

%% Data length resampler
function data1_interp = DataLengthResampler(data1, data2)
% Resample data1 to match the length of data2 using normalized index scaling

% Create normalized index for each dataset
norm1 = linspace(0, 1, length(data1));
norm2 = linspace(0, 1, length(data2));

% Interpolate data1 onto the normalized index of data2
% data1_interp = interp1(norm1, data1, norm2, 'makima');
if length(data1)<length(data2)
    data1_interp = interp1(norm1, data1, norm2, 'spline');
else
    data1_interp = interp1(norm1, data1, norm2, 'makima');
end
end
