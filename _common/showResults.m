function [] = showResults(targetTracks, targetExtents, estimatedTracks, estimatedExtents, measurementsCell, axisValues, visualizationMode)
numTargets = size(targetTracks,3);
numSteps = size(measurementsCell,1);
mycolors = [0.66,0.00,0.00; 0.00,0.30,0.70; 0.50,0.00,0.00];

for step = 1:numSteps
    
    if(visualizationMode == 0)
        tmp = numSteps;
    else
        tmp = step;
    end
    
    % plot the true target trajectory
    fig1 = figure(1);
    set(fig1, 'Position', [50, 50, 800, 700]);
    scatter([10000;10000],[10000;10000],50,mycolors(1,:),'x','LineWidth', 1.5)
    hold on
    for target = 1:numTargets
        scatter(targetTracks(1,tmp,target),targetTracks(2,tmp,target),50,mycolors(1,:),'x','LineWidth', 2)
        p = plot(targetTracks(1,1:tmp,target), targetTracks(2,1:tmp,target),'LineWidth', 1.5,'MarkerSize',25);
        set(p,'color',mycolors(1,:));
        if(~any(isnan(targetExtents(:,:,tmp,target)^2)))
            p = errorEllipse(targetExtents(:,:,tmp,target)^2,targetTracks(1:2,tmp,target));
            set(p,'color',mycolors(1,:));
        end
    end
    
    % plot all measurements
    measurements = measurementsCell{tmp};
    scatter(measurements(1,:),measurements(2,:),40,mycolors(3,:),'.','LineWidth', 1)

    
    numTracks = size(estimatedTracks,3);
    % plot the estimated target trajectory
    for target = 1:numTracks
        scatter(estimatedTracks(1,tmp,target),estimatedTracks(2,tmp,target),50,mycolors(2,:),'x','LineWidth', 2)
        p = plot(estimatedTracks(1,1:tmp,target), estimatedTracks(2,1:tmp,target),'LineWidth', 1.5,'MarkerSize',25);
        set(p,'color',mycolors(2,:));
        if(~any(isnan(estimatedExtents(:,:,tmp,target)^2)))
            p = errorEllipse(estimatedExtents(:,:,tmp,target)^2,estimatedTracks(1:2,tmp,target));
            set(p,'color',mycolors(2,:));
        end
        
        
        
    end
    hold off
    
    % adjust axis
    axis (axisValues);
    xlabel('x-coordinate [m]'), ylabel('y-coordinate [m]');
    pbaspect('manual')
    
    
    % pause the video if applicable
    if(step == 1 && visualizationMode ~= 0)
        pause
    end
    
    if(visualizationMode == 0)
        break
    elseif(visualizationMode >= 2)
        pause
    else
        pause(0.0001)
    end
    
end

end