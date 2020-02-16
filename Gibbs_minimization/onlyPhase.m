% This function returns a form of some original matrix (old) with only species in it of the specified phase (phaseNumber)
function new = onlyPhase(phaseNumber, old)           %Takes a given matrix and removes any values from that matrix that are not associated with.
    global l                                         % the specificed phase.
    new = [];                                        % The new version of the matrix with only the given phase.
    for temp = 1:length(old)                         % For each species,
        if l(temp)==phaseNumber                      %  if the species phase is the correct phase
            new = [new;old(temp)];                   %  add it to the new matrix (thus, the new matrix will only contain species of the
        end                                          %  correct phase.
    end
end