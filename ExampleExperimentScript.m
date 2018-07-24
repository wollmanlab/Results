%% Init
clear classes
close all
clc

%% load Experiment object
projectpth = '/data2/Images/anna/EGF_response_Curve'; 
Name = 'TestingErkPositiveFeedbackWithMMPinhibition'; 
E = Experiment(projectpth,Name); 

%% Create the relevant results objects
R(1) = MultiPositionSingleCellResults('/data2/Images/anna/EGF_response_Curve/EKARev_EGF_GM6001inhibition4_plus_FusarisetinA_2014May20'); 
R(2) = MultiPositionSingleCellResults('/data2/Images/anna/EGF_response_Curve/EKARev_EGF_GM6001inhibition6_plus_FusarisetinA_2014Jun09'); 


%% get figures and update their names to be unique
DoseResponseMatrices052014 = R(1).getFigure('DoseResponseMatrices'); 
DoseResponseMatrices052014.name = 'DoseResponseMatrices052014'; 

DoseResponseMatrices060914 = R(2).getFigure('DoseResponseMatrices');
DoseResponseMatrices060914.name = 'DoseResponseMatrices060914';

%% Add figures to the Experiment (they all must have unique names...)
E.addFigure(DoseResponseMatrices052014); 
E.addFigure(DoseResponseMatrices060914); 

%% create a new figure that compare the two datasets
ComparisonFigrue = Figure_(E,'compare_060914_and_052014'); 

%% use the resutls objects to get the data
Erk1 = R(1).getTimeseriesData('Erk'); 
Erk2 = R(1).getTimeseriesData('Erk'); 

%% plot new figure
plot([mean(Erk1) mean(Erk2)])

%% save figures
ComparisonFigrue.saveFigure; 

%% save Experiment
E.saveExperiment; 

