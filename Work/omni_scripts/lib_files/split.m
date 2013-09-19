function fields = split(string,delimiter)
%SPLIT Splits a string Perl style
%
% fields = split(string,delimiter)
%
% Default delimiter is '\W' (whitespace)
% Results are returned in the cell array fields 
%
% 07/14/04 Markus Herrgard

if (nargin < 2)
    delimiter = '\W';
end

[start_ind,end_ind] = regexp(string,delimiter);

if (~isempty(start_ind))
cnt = 0;
for i = 1:length(start_ind)+1
    if (i == 1)
        if (end_ind(i) > 1)
            cnt = cnt + 1;
            fields{cnt} = string(1:end_ind(i)-1);    
        end
    elseif (i == length(start_ind)+1)
        if (start_ind(i-1) < length(string))
            cnt = cnt + 1;
            fields{cnt} = string(start_ind(i-1)+1:end);
        end
    else
        cnt = cnt + 1;
        fields{cnt} = string(start_ind(i-1)+1:end_ind(i)-1);
    end
end
else
    fields{1} = string;
end

fields_out = {};
cnt = 0;
for i = 1:length(fields)
    if (~isempty(fields{i}))
        cnt = cnt+1;
        fields_out{cnt} = fields{i};
    end
end
fields = fields_out;