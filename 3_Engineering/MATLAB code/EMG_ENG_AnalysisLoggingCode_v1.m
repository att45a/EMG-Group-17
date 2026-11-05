s = serialport("COM6", 115200);  % change: type serialportlist to find it
pause(2); flush(s);               % let the Arduino restart, clear old data
readline(s);                      % throw away the first (possibly half) line

fs = 1000;                        % must match the Arduino (1000 Hz)
N = 10 * fs;                      % 30 seconds
Signal = zeros(1, N);
disp('Recording')
for k = 1:N
    Signal(k) = str2double(readline(s));
end
clear s
disp('Done')

Time = (1:N) / fs;
save('op_test2.mat', 'Time', 'Signal');   % change the name for each take
plot(Time, Signal)