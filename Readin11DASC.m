%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
% Auto-generated by MATLAB on 30-Jan-2023 00:11:48


%% Set source location

sourcefolder = ["/Users/marthanarinemrihavenith/Desktop/Manuscripts/Submitted/2023-10 Breathwork CO2/Source Excels/11DASC.xlsx"];


%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:K62";

% Specify column names and types
opts.VariableNames = ["ExpofUnityEoU", "SpiritualExpSE", "BlissfulStateBS", "InsightfullnessI", "DisembodimntD", "ImpairedControlandCognitionICaC", "AnxietyA", "CompleyImageryCI", "ElementaryImageryEI", "AudioVisualSynesthesiaeAVS", "ChangedMeanungofPerceptsCMoP"];
opts.SelectedVariableNames = ["ExpofUnityEoU", "SpiritualExpSE", "BlissfulStateBS", "InsightfullnessI", "DisembodimntD", "ImpairedControlandCognitionICaC", "AnxietyA", "CompleyImageryCI", "ElementaryImageryEI", "AudioVisualSynesthesiaeAVS", "ChangedMeanungofPerceptsCMoP"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
tbl = readtable(sourcefolder, opts, "UseExcel", false);

%% Convert to output type
DASC_ExperienceOfUnity = tbl.ExpofUnityEoU;
DASC_SpiritualExperience = tbl.SpiritualExpSE;
DASC_BlissfulState = tbl.BlissfulStateBS;
DASC_Insightfulness = tbl.InsightfullnessI;
DASC_Disembodiment = tbl.DisembodimntD;
DASC_ImpairedControlandCognition = tbl.ImpairedControlandCognitionICaC;
DASC_Anxiety = tbl.AnxietyA;
DASC_ComplexImagery = tbl.CompleyImageryCI;
DASC_ElementaryImagery = tbl.ElementaryImageryEI;
DASC_AudioVisualSynesthesia = tbl.AudioVisualSynesthesiaeAVS;
DASC_ChangedMeaningofPercepts = tbl.ChangedMeanungofPerceptsCMoP;

DASCmatrix =[DASC_ExperienceOfUnity, DASC_SpiritualExperience, DASC_BlissfulState,DASC_Insightfulness,DASC_Disembodiment,DASC_ImpairedControlandCognition,DASC_Anxiety, DASC_ComplexImagery, DASC_ElementaryImagery, DASC_AudioVisualSynesthesia, DASC_ChangedMeaningofPercepts];
DASCnames = {'Unity';'Spiritual';'Bliss';'Insight';'Disembodiment';'Impaired Control'; 'Anxiety';'Complex Imagery'; 'Simple Imagery'; 'Synaesthesia'; 'Changed meaning'};
%% Clear temporary variables
clear opts tbl DASC_* 