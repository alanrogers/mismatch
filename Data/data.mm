
% Results will include the reference sequence (line 1).
% Expecting 118 subjects X 312 sites = 36816 sites in all

% 118 subjects 312 sites
% mean diff = 6.105751.  segregating sites = 24

% A. FORMATTED FOR MMEST:
sampsize = 118 ; % number of subjects
segregating = 24 ; % number of segregating sites
nsites = 312 ; % number of sites
%mismatch distribution:
histogram = 1271 1141 1265 773 128 0 0 0 0 0 0 0 10 62 297 1222 361 258 115 ;

% B. NORMALIZED HISTOGRAM:
% histogram = 0.184123 0.165290 0.183254 0.111980
%   0.018543 0.000000 0.000000 0.000000 0.000000
%   0.000000 0.000000 0.000000 0.001449 0.008982
%   0.043025 0.177024 0.052296 0.037375 0.016659

%C. NORMALIZED HISTOGRAM IN PicTeX FORMAT:
%\plot
%   0 0.184123 1 0.165290 2 0.183254 3 0.111980
%   4 0.018543 5 0.000000 6 0.000000 7 0.000000 8 0.000000
%   9 0.000000 10 0.000000 11 0.000000 12 0.001449 13 0.008982
%   14 0.043025 15 0.177024 16 0.052296 17 0.037375 18 0.016659
%/
