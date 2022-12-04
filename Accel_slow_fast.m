function acceleration=Accel_slow_fast(time_start, time_span, currentTime, distance_change)

mu=time_start+time_span/2;
sigma=time_span/8;
acceleration=-distance_change*1/(sqrt(2*pi*sigma^2))*exp(-(currentTime-mu)^2/2/sigma^2)*(-currentTime+mu)/sigma^2;