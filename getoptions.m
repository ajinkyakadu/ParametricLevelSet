function v = getoptions(options, name, v, mandatory)

% getoptions - retrieve options parameter



if nargin<4
    mandatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mandatory
    error(['You have to provide options.' name '.']);
end 