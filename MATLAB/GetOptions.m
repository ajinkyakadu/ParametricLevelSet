function optionValue = GetOptions(optionsStruct, optionName, defaultValue, isMandatory)
% GETOPTIONS - Retrieve specified option from options struct
%
% Syntax:
%   optionValue = GetOptions(optionsStruct, optionName, defaultValue, isMandatory)
%
% Inputs:
%   optionsStruct - A structure containing options and their respective values.
%   optionName - The name of the option to retrieve.
%   defaultValue - The default value of the option, used if the option is not present in the options struct.
%   isMandatory - Boolean indicating if the option is mandatory (default: false).
%
% Outputs:
%   optionValue - The retrieved option value.
%
% Author: 
%   Ajinkya Kadu
%
% Date: 
%   September 2016

if nargin < 4
    isMandatory = false;
end

if isfield(optionsStruct, optionName)
    optionValue = eval(['optionsStruct.' optionName ';']);
elseif isMandatory
    error(['Mandatory option ' optionName ' not found in options struct.']);
else
    optionValue = defaultValue;
end
