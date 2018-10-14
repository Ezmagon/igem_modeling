%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% increaseVersion(bumpType)
%
% Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function increaseVersion(bumpType)

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if ~strcmp(currentBranch,'master')
    error('ERROR: not in master')
end

%Bump version number:
oldModel   = load('../ModelFiles/mat/yeastGEM.mat');
oldVersion = oldModel.model.description;
oldVersion = oldVersion(strfind(oldVersion,'_v')+2:end);
oldVersion = str2double(strsplit(oldVersion,'.'));
newVersion = oldVersion;
switch bumpType
    case 'major'
        newVersion(1) = newVersion(1) + 1;
        newVersion(2) = 0;
        newVersion(3) = 0;
    case 'minor'
        newVersion(2) = newVersion(2) + 1;
        newVersion(3) = 0;
    case 'patch'
        newVersion(3) = newVersion(3) + 1;
    otherwise
        error('ERROR: invalid input. Use "major", "minor" or "patch"')
end
newVersion = num2str(newVersion,'%d.%d.%d');

%Check if history has been updated:
fid     = fopen('../history.md','r');
history = fscanf(fid,'%s');
fclose(fid);
if ~contains(history,['yeast' newVersion ':'])
    error('ERROR: update history.md first')
end

%Load model:
initCobraToolbox
model = readCbModel('../ModelFiles/xml/yeastGEM.xml');

%Include tag and save model:
model.modelID = ['yeastGEM_v' newVersion];
saveYeastModel(model,false)

%Allow .mat & .xls storage:
copyfile('../.gitignore','backup')
fin  = fopen('backup','r');
fout = fopen('../.gitignore','w');
still_reading = true;
while still_reading
  inline = fgets(fin);
  if ~ischar(inline)
      still_reading = false;
  elseif ~startsWith(inline,'*.mat') && ~startsWith(inline,'*.xlsx')
      fwrite(fout,inline);
  end
end
fclose('all');
delete('backup');

%Store model as .mat:
save('../ModelFiles/mat/yeastGEM.mat','model');

%Convert to RAVEN format and store model as .xlsx:
model = ravenCobraWrapper(model);
exportToExcelFormat(model,'../ModelFiles/xlsx/yeastGEM.xlsx');

%Update version file:
fid = fopen('../version.txt','wt');
fprintf(fid,newVersion);
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%