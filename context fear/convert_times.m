timeStr = '43:49'; % Example MM:SS string
minutesSeconds = sscanf(timeStr, '%d:%d'); % Extract minutes and seconds
totalSeconds = minutesSeconds(1) * 60 + minutesSeconds(2); % Convert to seconds
disp(totalSeconds);
%%

timeStr = '12:00'; % Example MM:SS string
minutesSeconds = sscanf(timeStr, '%d:%d'); % Extract minutes and seconds
totalSeconds = minutesSeconds(1) * 60 + minutesSeconds(2); % Convert to seconds
disp(totalSeconds);

%%
timeStr = '14:00'; % Example MM:SS string
minutesSeconds = sscanf(timeStr, '%d:%d'); % Extract minutes and seconds
totalSeconds = minutesSeconds(1) * 60 + minutesSeconds(2); % Convert to seconds
disp(totalSeconds);

840/2629.533