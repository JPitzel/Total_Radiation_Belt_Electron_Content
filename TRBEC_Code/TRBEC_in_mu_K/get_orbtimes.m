
%computes half orbits given position data
%NOTE: marks times independently of units for position (some of the MagEIS data is given in R_E and others in km)
function orbtimes = get_orbtimes(filenames)
orbtimes = [];
for file_idx = 1:length(filenames)
    [data_orbtimes,~] = spdfcdfread(filenames{file_idx},'Variables',{'Epoch','Position'});
    epoch = data_orbtimes{1};
    position = data_orbtimes{2};
    distance = sqrt(sum(position.^2,2));
    
    %remove data that is too large
    epoch(distance>6*6371) = [];
    distance(distance>6*6371) = [];
    
    %if average is less than 10 R_E, convert to units of km
    %this fixes an issue where some versions of the MagEIS data (e.g. v8.4.0) use units of R_E instead of km as expected
    if mean(distance) < 10
        distance = distance*6371;
    end
    
    %smooth to avoid discontinuities that occur at apogee and perigee (i.e. there were instances where expected local min appeared as two local mins and one local max)
    distance = smooth(distance,21);
    
    %sort out local min and max that indicate satellite at perigee and apogee
    local_max_idx = islocalmax(distance);
    local_min_idx = islocalmin(distance);
    orbtimes_i = [epoch(local_max_idx);epoch(local_min_idx)];%combine and sort by occurance
    orbtimes_i = sortrows(orbtimes_i);

    orbtimes = cat(1,orbtimes,orbtimes_i);

    %plot distance vs time for troubleshooting
    %plot(epoch,distance)
    %hold on
end

%once the orbtimes are determined, remove orbits that are too short (some still occur even with smoothing done above)
orbperiod = diff(orbtimes);
idxs = find(orbperiod<0.95*mean(orbperiod));%find indices where difference between orbtimes is less than 95 percent of mean half orbit period
orbtimes(idxs+1) = [];%remove the identified orbtimes that were too short
end
