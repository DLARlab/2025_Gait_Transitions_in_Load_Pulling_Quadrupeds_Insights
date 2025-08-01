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
