function neuron = t2n_setionconcentration(neuron,str)
% This function modifies the Nernst potential or the extra- and
% intracellular ion concentrations (which then leads to different Nernst
% potential) according to specific experiments in literature. The list can
% be extended.
%
% INPUTS
% neuron    neuron structure containing ion channel densities 
% str       string defining the experimental conditions from literature.
%           Current cases: 'Krueppel' from Krueppel et al. (2011) Neuron
%                          'Mongiat'&'Mongiat2' from Mongiat et al (2009) PlOS One
%                          'Riazanski'&'SH07' from Schmidt-Hieber et al (2007) JNS; and Riazanski et al (2001) J. of Physiology
%                          'SH08' from Schmidt-Hieber et al (2008) Journal of Physiology
%                          'Stocca' from Stocca et al (2008) Journal of Physiology, (used for the adult mice)
%
% OUTPUTS
% neuron    updated neuron structure
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

switch str
    case 'Krueppel'
        % from Krueppel et al. (2011) Neuron
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-90);
            neuron.mech{t}.all.na_ion = struct('ena',60);
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);
        end
    case 'Mongiat'
        % calculated from concentrations in Mongiat et al (2009) PlOS One
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-93);               
            neuron.mech{t}.all.na_ion = struct('ena',87.76);  
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048); 
        end
    case 'Mongiat2'
        % calculated from concentrations in Mongiat et al (2009) PlOS One
        % use this version if nakbuffer is available too.
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ko0',4,'ki0',140);
            neuron.mech{t}.all.na_ion = struct('nao0',156,'nai0',5);
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);
        end
    case {'Riazanski','SH07'}
        % Schmidt-Hieber et al (2007) JNS; and Riazanski et al (2001) J. of
        % Physiology
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-103.2);
            neuron.mech{t}.all.na_ion = struct('ena',82.59);   
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);   
        end
    case 'SH08' 
        % Schmidt-Hieber et al (2008) Journal of Physiology
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-103.2);
            neuron.mech{t}.all.na_ion = struct('ena',81.07);   
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);   
        end
    case 'Stocca'
        % Stocca et al (2008) Journal of Physiology, for the adult animals
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-103.2);
            neuron.mech{t}.all.na_ion = struct('ena',78.78);   
            neuron.mech{t}.all.ca_ion = struct('cao0',0.5,'cai0',0.000048);   
        end
    otherwise
        error('String "%s" not found in the list',str)
end