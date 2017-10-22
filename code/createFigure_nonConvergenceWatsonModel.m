function createFigure_nonConvergenceWatsonModel( varargin )
% createFigure_nonConvergenceWatsonModel( varargin )
%
% This routine examines the performance of the Watson 2014 equations
% when applied to the Drasdo 2007 model for calculating the radial
% displacement of retinal ganglion cells.
%
% OPTIONS
%   sampleResolutionDegrees - The calculations are performed across a
%       regular sampling of eccentricity. This param defines sample
%       resolution. The sample resolution must be sufficient fine so that
%       the cumulative is an accurate estimate of the integral.
%   maxModeledEccentricity - The eccentricity extent of the model.
%   meridianAngleResolutionDeg - The resolution across polar angle for
%       which displacements are calculated.
%   verbose - Do we give you the text?
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',90,@isnumeric);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('pathToPlotOutputDir','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})


%% Setup
% Prepare the regular eccentricity support base
regularSupportPosDeg = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% Prepare the set of meridian angles for which we will calculate
% displacement
meridianAngles = 0:p.Results.meridianAngleResolutionDeg:(360-p.Results.meridianAngleResolutionDeg);

% Prepare a figure
figHandle = figure();

%% Loop over the meridians
for mm = 1:length(meridianAngles)
    
    %% mRF_cumulative function
    % We build a function that returns the cumulative mRF density.
    
    % Start with Watson (2014) eq 8.
    mRFDensityOverRegularSupport = calcWatsonMidgetRFDensityByEccen(regularSupportPosDeg, meridianAngles(mm));
    % Define the cumulative sum of mRF density
    mRF_cumulative = calcCumulative(regularSupportPosDeg, mRFDensityOverRegularSupport);
    
    
    %% mRGC_cumulative function
    % We build a function that returns the cumulative mRGC density
    % Obtain a spline fit to the empirical RGC density data of Curcio 1990
    RGCDensityFit = getSplineFitToRGCDensity(meridianAngles(mm));
    % Create a function that returns mRGC density as a function of
    % RGC density, with the transform defined by Watson (2014) Eq 7 (which
    % itself is taken from Drasdo 2017)
    mRGCDensityOverRegularSupport = ...
        RGCDensityFit(regularSupportPosDeg)' .* ...
        calcDrasdoMidgetFractionByEccen(regularSupportPosDeg,0.8928,41.03);
    % Define the cumulative sum of mRGC density
    mRGC_cumulative = calcCumulative(regularSupportPosDeg, mRGCDensityOverRegularSupport);
    
    % Calculate and store displacement and cumulative functions
    mRGC_cumulativeEachMeridian(mm,:)=mRGC_cumulative;
    mRF_cumulativeEachMeridian(mm,:)=mRF_cumulative;
    rgcDisplacementEachMeridian(mm,:)=calcDisplacement(regularSupportPosDeg, mRGC_cumulative, mRF_cumulative);
    
    % Report the results for this meridian
    if p.Results.verbose
        zeroPoints = find(rgcDisplacementEachMeridian(mm,:)==0);
        convergenceIdx = find(regularSupportPosDeg(zeroPoints) > 2,1);
        if isempty(convergenceIdx)
            outLine = ['Polar angle: ' num2str(meridianAngles(mm)) ', max RGC displacement: ' num2str(max(rgcDisplacementEachMeridian(mm,:))) ', found convergence: FAILED TO CONVERGE\n'];
            fprintf(outLine);
        else
            convergenceEccen(mm) = regularSupportPosDeg(zeroPoints(convergenceIdx));
            outLine = ['Polar angle: ' num2str(meridianAngles(mm)) ', max RGC displacement: ' num2str(max(rgcDisplacementEachMeridian(mm,:))) ', found convergence: ' num2str(convergenceEccen(mm)) '\n'];
            fprintf(outLine);
        end
    end
    
    % plot the displacement
    subplot(length(meridianAngles),2,mm*2);
    plot(regularSupportPosDeg,rgcDisplacementEachMeridian(mm,:),'-r')
    axis off;
    ylim([-.5 3.0]);
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('RGC displacement [deg]');
    end
    
    % Plot the cumulative functions
    subplot(length(meridianAngles),2,mm*2-1);
    plot(regularSupportPosDeg,mRGC_cumulativeEachMeridian(mm,:),'-k')
    axis off;
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('cells per sector');
    end
    hold on
    plot(regularSupportPosDeg,mRF_cumulativeEachMeridian(mm,:),'-b')
    ylim([0 8e5]);
    hold off
    drawnow
    
end % loop over meridians

if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'TestWatsonConvergence.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


end % createFigure_nonConvergenceWatsonModel


%% LOCAL FUNCTIONS

function displaceInDeg = calcDisplacement(regularSupportPosDeg, countPerRingRGC, countPerRingRF)

% Determine the sample resolution by a difference operation
tmp = diff(regularSupportPosDeg);
sampleResolutionDegrees = tmp(1);
% Measure the displacement (in degrees). First, for each cumulative RGC
% density value, identify the array index of the first
% cumulative RF density value that is equal or greater.
displaceInSamples=arrayfun(@(x) find(countPerRingRF>=x,1), countPerRingRGC,'UniformOutput',false);
% Now some array operations to get these values out of cells and in to a
% numeric vector
emptyCells = find(cellfun(@(x) isempty(x), displaceInSamples));
displaceInSamples(emptyCells)={NaN};
displaceInSamples=cell2mat(displaceInSamples(cellfun(@(x) ~isempty(x), displaceInSamples)));
% The displacement in array samples is the array index of each cumulative
% RGC density measurement, minus the array index of the matching RF density
% measurement. This difference is then multiplied by the sample resolution
% to get the displacement in degrees.
displaceInDeg = ((1:1:length(displaceInSamples))-displaceInSamples ) * sampleResolutionDegrees;
% Zero out negative values after the point of convergence
displaceInDeg(find(displaceInDeg < 0,1):end)=0;

end

