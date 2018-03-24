%% 440305585
% AERO4701
%
% Calculates the date and time given an Epoch year and decimal day
%
% inputs:   sat, Epoch year, Epoch decimal days
% outputs:  [month, day, hour, min, sec]

function sat = date_time_TLE(sat, Epoch_year, Epoch_days)

    %%  year
    
    if (Epoch_year < 57)          % first satellite appeared in 1957
        sat.year = 2000 + Epoch_year; 
    else
        sat.year = 1900 + Epoch_year; 
    end
    
    %% month and day
        
    if (leapyear(sat.year))
        days_in_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    else
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    end
    
    month = 1;
    while(Epoch_days > days_in_month(month))
        Epoch_days = Epoch_days - days_in_month(month);
        month = month + 1;
    end
    sat.month = month;
    sat.day = floor(Epoch_days);
    Epoch_hours = (Epoch_days - sat.day)*24;
    
    %% time
    
    sat.hour = floor(Epoch_hours);
    Epoch_mins = (Epoch_hours - sat.hour)*60;
    sat.min = floor(Epoch_mins);
    sat.sec = floor((Epoch_mins - sat.min)*60);
    
end