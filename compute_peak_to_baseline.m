function [R,me,ml,ie,il] = compute_peak_to_baseline(x, y)
% Computes the peak to baseline ration according to Phoka et al.2010
% [R,me,ml,ie,il] = compute_peak_to_baseline(x, y)

[me,ie] = find_max_local_extrema(y(x<0.5 & x>=0));
if ~isempty(me)
    ie = ie + find(x>=0,1,'first');
else
    fprintf(1,'Not enough points or data with no peaks...\n')
    fprintf(1,'Using absolute maxima!\n')
    [me,ie] = max(y(x<0.5 & x>=0));
end
[ml,il] = find_max_local_extrema(y(x>=0.5&x<1));

if ~isempty(ml)
    il = il + find(x>=0.5,1,'first');
else
    fprintf(1,'Not enough points or data with no peaks...\n')
    fprintf(1,'Using absolute maxima!\n')
    [ml,il] = max(y(x>=0.5&x<1));
    il = il + find(x>=0.5,1,'first');

end

if ~isempty(me) & ~isempty(ml)
    R = abs(me-ml)./(abs(me)+abs(ml));
else
    fprintf(1,'Check data...')
    keyboard
    
end



function [m, ii] = find_max_local_extrema(y)

[m, idx] = findpeaks(abs(y));
[m,ii] = max(m);
ii = idx(ii);
m = y(ii);   
% fig = figure(10);
% plot(y),hold all
% plot(abs(y))
% plot(ii,m,'o')
% pause(1)
% close(fig)
