function add_stack(name)
if (ispc == 1)
    filename = strcat('../Models/', name, '/Code/NEURON/TMS_script.ps1');
    filename2 = strcat('../Models/', name, '/Code/NEURON/save_locations.ps1');
    [~,outp] = system('where nrniv');
    if isempty(outp) || ~isempty(regexp(outp,'not found','ONCE'))  || ~isempty(regexp(outp,'no nrniv','ONCE')) || ~isempty(regexp(outp,'not find','ONCE')) || strcmp(outp,sprintf('\n')) || isempty(regexp(outp,'nrniv','ONCE'))
        t2npath = fileparts(which('t2n.m'));
        fid = fopen(fullfile(t2npath,'nrniv_win.txt'),'r');
        nrnivPath = fread(fid,'*char')';
        fclose(fid);
    else
        nrnivPath = 'nrniv';
    end
else
    filename = strcat('../Models/', name, '/Code/NEURON/TMS_script.sh');
    filename2 = strcat('../Models/', name, '/Code/NEURON/save_locations.sh');
    nrnivPath = 'nrniv';
end

fileID = fopen(filename,'w');
fprintf(fileID,[nrnivPath ' -NSTACK 100000 -NFRAME 20000 -Py_NoSiteFlag TMS_script.hoc']);
fclose(fileID);
fileID = fopen(filename2,'w');
fprintf(fileID,[nrnivPath ' -NSTACK 100000 -NFRAME 20000 -Py_NoSiteFlag save_locations.hoc']);
fclose(fileID);