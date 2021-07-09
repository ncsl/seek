function exportTSV(elec, fname)
    % Creates a `.tsv` file from Fieldtrip `cfg.elec` data structure.
    %
    % Meant for saving data in `*electrodes.tsv` BIDS format.
    %
    % Inputs:
    %   elec    : a struct with `label`, and `chanpos` fields.
    %   file    : The filepath to save.

    if ~contains(fname, 'electrodes.tsv')
        error(strcat('Filename should be of the form "*electrodes.tsv". ',... 
            'The filename does not have electrodes.tsv substring inside.'));
    end
    
    % create separator
    sep = regexp(elec.label,'\d');
    
    % get the vector of x, y, and z coordinates
    xcoord = elec.chanpos(:,1);
    ycoord = elec.chanpos(:,2);
    zcoord = elec.chanpos(:,3);

    % write file to 'name'
    fname
    fid = fopen(fname, 'w');
    fprintf(fid,'name\tx\ty\tz\n');
    for i = 1:length(elec.label)
        elecName = strcat(elec.label{i}(1:sep{i}-1),"'",elec.label{i}((sep{i}:length(elec.label{i}))));
        fprintf(fid, '%s\t%-.3f\t%-.3f\t%-.3f\n', elecName, xcoord(i), ycoord(i), zcoord(i));
    end
    fclose(fid)
end
