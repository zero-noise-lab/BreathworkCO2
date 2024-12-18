%% Import data from MEQ spreadsheet
%
% Auto-generated by MATLAB on 31-May-2022 23:27:13

%% Set source location

sourcefolder = ["/Users/marthanarinemrihavenith/Desktop/Manuscripts/Submitted/2023-10 Breathwork CO2/Source Excels/MEQ.xlsx"];

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:G34";

% Specify column names and types
opts.VariableNames = ["TranscendenceP", "PositiveMoodP", "IneffabilityI", "Mystical", "TOTAL","HolConMEQ","ActPasMEQ"];
opts.SelectedVariableNames = ["TranscendenceP", "PositiveMoodP", "IneffabilityI", "Mystical", "TOTAL","HolConMEQ","ActPasMEQ"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Import the data
MEQ = readtable(sourcefolder, opts, "UseExcel", false);

% Organize data
holconSURVEYvector = MEQ.HolConMEQ;
breathSURVEYvector = MEQ.ActPasMEQ+1;
MEQmatrix = [MEQ.TranscendenceP, MEQ.PositiveMoodP, MEQ.IneffabilityI,MEQ.Mystical,MEQ.TOTAL];
MEQnames = {'Transcendence'; 'Positive Mood'; 'Ineffability'; 'Mysical Experience'; 'Total score'};

%% Clear temporary variables
clear opts MEQ