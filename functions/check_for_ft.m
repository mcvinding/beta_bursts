function check_for_ft()
% Check for FieldTrip in paths

if exist('ft_defaults','file')
    ft_defaults
else
   fprintf('Could not find ft_defaults in path\nPleas add FieldTrip folder to path and try again!')
end

%END
   