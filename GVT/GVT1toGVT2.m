% Title: Convert GVT1 Data to GVT2 Format
% Author: Anthony Su
% Date: 2023-07-14

function new_obj = GVT1toGVT2(old_obj_filename)

% load old object
old_obj = load(old_obj_filename).H;

% import rate and sensors
new_obj.rate = old_obj(1).rate;
new_obj.sensors = old_obj(1).sensors;

% correct hammer (sensor 1) location
temp = split(string(old_obj_filename),'impact_'); % cut off before hammerLoc
temp = split(temp(2),'_sens'); % cut off after hammerLoc
hammerLoc = temp(1);
new_obj.sensors(1).Location = {char(hammerLoc)};


% import data

for idx_impact = 1:length(old_obj)
    H = [];
    H.data = old_obj(idx_impact).data;

    % find impact
    for idx = 1:width(H.data)
        [pvs_,pis_] = findpeaks(H.data{:,idx});
        [pv_,pi_] = max(pvs_);
        peakVals(idx) = pv_;
        peakIdxs(idx) = pis_(pi_);
    end
    H.peakIdxs = peakIdxs;

    % assign to new_obj
    new_obj.H(idx_impact) = H;
end

end