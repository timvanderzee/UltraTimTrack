function parms = adjust_hough_parameters(parms)

%% Init variables
adj_parms = [parms.fas.range parms.apo.super.maxangle parms.apo.deep.maxangle]; %get the 4 parms I can change
%% initialize text var for the inputs
% Create a figure window
fig = figure('Units', 'normalized', 'Position', [0.4 0.4 0.3 0.5], ...
    'Name', 'Set TimTrack settings', 'NumberTitle', 'off');

text_prompt = {'Min. fascicle angle (deg):', 'Max. fascicle angle (deg):', ...
    'Max. superficial aponeurosis angle (deg):', 'Max. deep aponeurosis angle (deg):',...
    'The updated parameters will be saved, reloaded, and automatically applied in future UltraTimTrack sessions'};
fields = {'fas_min', 'fas_max', 'sup_apo', 'deep_apo','info'};
num_fields = numel(fields);
inputs = gobjects(num_fields, 1);

%% Create input fields
for i = 1:num_fields-1
    uicontrol('Style', 'text', 'String', text_prompt{i}, ...
        'Units', 'normalized', 'Position', [0.1, 0.85 - (i-1) * 0.18, 0.5, 0.08], ...
        'FontSize', 14, 'HorizontalAlignment', 'left');

    inputs(i) = uicontrol('Style', 'edit', 'String', num2str(adj_parms(i)), ...
        'Units', 'normalized', 'Position', [0.6, 0.85 - (i-1) * 0.18, 0.3, 0.08], ...
        'FontSize', 14, 'Callback', @(src, ~) update_param(i, src));
end

% Info Text (aligned properly under inputs)
uicontrol('Style', 'text', 'String', text_prompt{end}, ...
    'Units', 'normalized', 'Position', [0.1, 0.20, 0.8, 0.08], ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Save Button (Centered)
uicontrol('Style', 'pushbutton', 'String', 'Save', ...
    'Units', 'normalized', 'Position', [0.35, 0.05, 0.3, 0.1], ...
    'FontSize', 14, 'FontWeight', 'bold', 'Callback', @okButtonCallback);

%% Callback function to update parameter values
    function update_param(field, src)
        value = str2double(get(src, 'String'));
        if isnan(value) || isempty(value) || (value >= 90)
            set(src, 'String', num2str(adj_parms(field)));
            errordlg('Invalid input. Please enter a numeric value <90.', 'Error', 'modal');
        else        
            adj_parms(field) = value;
        end
    end

%% OK Button Callback and update parms to the struct before returning it
    function okButtonCallback(~, ~)
        parms.fas.range = adj_parms(1:2);
        parms.apo.super.maxangle = adj_parms(3);
        parms.apo.deep.maxangle = adj_parms(4);
        close(fig);
    end

% Wait for the figure to close before returning updated parameters
uiwait(fig);

end
