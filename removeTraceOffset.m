function [x, offset] = removeTraceOffset(x, idx)
% [x] = removeTraceOffset(x, idx)
%
% Removes the offset of a trace.
% x is the trace. This function acts on the rows of x.
% idx [A,B] is an array containing the initial and end index used on the
% averaging.
% Offset is the offset used to compute
%
iremove_offset = @(y)(y - mean(y(idx(1):idx(2))));

if iscolumn(x) || isrow(x)
    offset = mean(y(idx(1):idx(2)));
    x = iremove_offset(x);
else
    m = size(x,1);
    offset = zeros(m,1);
    for k = 1:m
        offset(k) = mean(x(k,idx(1):idx(2)));
        x(k,:) = iremove_offset(x(k,:));
    end
end
