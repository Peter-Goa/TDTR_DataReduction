function [ max_amp ] = GetAmp(td,amp) 
global norm_time
%%%%%%%     norm_timeのときのampの値を返す

for i = 1:length(td)
    if td(i) >= norm_time
        a = amp(i);
        break;
    end
end
max_amp = a;
end

