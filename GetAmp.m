function [ max_amp ] = GetAmp(td,amp) 
% the time at which the value will be used to normalize amplitude of data [s]
norm_time = 1.00E-10;
for i = 1:length(td)
    if td(i) >= norm_time
        a = amp(i);
        break;
    end
end
max_amp = a;
end

