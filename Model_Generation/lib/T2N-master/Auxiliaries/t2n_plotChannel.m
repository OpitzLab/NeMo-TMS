function t2n_plotChannel(neuron,mcondChannel,region,outputFolder,options)
% This function plots the activation and inactivation dynamics of an voltage- 
% dependent ion channel. The fits on the activation and inactivation time are
% mono-exponential, hence they should be handled with care in cases the
% channel has more complex kinetics.
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms
% mcondChannel     the name of the maximum conductance parameter of a
%                   channel, as it is in NEURON and NMODL (e.g. gbar_hh)
% region            name of the region in neuron_orig from which
%                   specifications for the channel should be taken from. If
%                   not provided, t2n takes the first specification it can
%                   find within a region
% outputFolder      (optional) folder where pdfs should be saved to
% options           string with more options that can be concatenated:
%                   -h          set if channel is hyperpolarization activated
%                   -si or -sa  set if channel inactivation or activation is very slow (dt and duration increased)
%                   -fi or -fa  set if channel inactivation or activation is very fast (dt and duration decreased)
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************
if ~exist('options','var') || isempty(options)
    options = '' ;
end
if ~exist('outputFolder','var') || isempty(outputFolder)
    outputFolder = '' ;
end
if ~exist(outputFolder,'dir')
    mkdir(outputFolder)
end
height = 4;
width = 4;
channel = strsplit(mcondChannel,'_'); mcond = channel{1};channel = channel{2};
cond_channel = strrep(mcondChannel,'bar','');

tree{1} = t2n_testComp;   % get the test compartment morphology

neuronTest.params = neuron.params;
neuronTest.params.nseg = 1;
neuronTest.params.cvode = 0;
neuronTest.mech{1}.all.(channel) = struct(mcond,0.1); % put the channel in all (so the one) compartments

if exist('region','var') && ~isempty(region)  % get the ion channel specification from the region
    neuronTest.mech{1}.all.(channel) = neuron.mech{1}.(region).(channel);
else   % else get it from the first region which incorporates the ion channel
    fields = fieldnames(neuron.mech{1});
    for f = 1:numel(fields)
        if isfield(neuron.mech{1}.(fields{f}),channel)
            neuronTest.mech{1}.all.(channel) = neuron.mech{1}.(fields{f}).(channel);
            region = fields{f};
            break
        end
    end
end

neuronTest.record{1}.cell = struct('node',1,'record',cond_channel);

neuronTest.params.cvode = 0;

if ~isempty(strfind(options,'-sa'))
    afac = double(cell2mat(textscan(options,'-sa%d')));
    if isempty(afac)
        afac = 1;
    end
    neuronTest.params.dt = 0.25*afac;
    dur = 50;
    actdur =  5000*afac;
elseif ~isempty(strfind(options,'-fa'))
    afac = double(cell2mat(textscan(options,'-fa%d')));
    if isempty(afac)
        afac = 1;
    end
    neuronTest.params.dt = 0.005/afac;
    dur = 10;
    actdur =  25/afac;
else
    neuronTest.params.dt = 0.025;
    dur = 50;
    actdur = 100;
end

% activation measurement
if ~isempty(strfind(options,'-h'))
    holding_voltage = 0;
    amp = -150:5:0;
else
    holding_voltage = -120;
    amp = -120:5:50;
end
[~,out] = t2n_voltSteps(neuronTest,tree,amp,[dur actdur dur],holding_voltage);
maxG = zeros(numel(amp),1);
tauA = maxG;
figure;hold all;
xlabel('Time [ms]')
ylabel('Conductance of channel')
for a = 1:numel(amp)
    thisVec = out{a}.record{1}.cell.(cond_channel){1};
    [maxG(a),indmaxG] = max(thisVec(out{a}.t>dur & out{a}.t<dur+actdur));
    
    indStart = find(out{a}.t > dur,1,'first');
    indmaxG = indmaxG + indStart -1;   % correct indmaxG
    if indmaxG - indStart > 20
        indmaxG = find(thisVec(indStart:indmaxG) > maxG(a)*0.99,1,'first')+indStart-1; % reduce length of fit vector for cases where there is nearly no slope towards maxG
    end
    % get part of vector which is linear, i.e. calculate derivative
    % and check where it deviates from the median by more than one
    % s.e.m.
    firstDiff = diff(log(-thisVec(indStart:indmaxG-1)+thisVec(indmaxG)),1,1);
    indStartNew = find(abs(firstDiff-median(firstDiff))<=(std(firstDiff)/sqrt(numel(firstDiff))*1.01),1,'first')+indStart;
    indEndNew = find(abs(firstDiff-median(firstDiff))<=(std(firstDiff)/sqrt(numel(firstDiff))*1.01),1,'last')+indStart+1;
    if indEndNew == indmaxG   % avoid zeros in log
        indEndNew = indEndNew-1;
    end
    if numel(indStartNew:indEndNew) > 1
        [ft,~] = polyfit(out{a}.t(indStartNew:indEndNew),log(-thisVec(indStartNew:indEndNew)+thisVec(indmaxG)),1);
        if imag(ft(1)) ~= 0
            ft = [NaN NaN];
        end
        plot(out{a}.t,thisVec)
        plot(out{a}.t(indStartNew:indEndNew),-exp(ft(1)*out{a}.t(indStartNew:indEndNew)+ft(2))+thisVec(indmaxG),'r--')
        tauA(a) = -1/ft(1);
        if abs(-1/ft(1)) > 1e6  % delete nonsense fits
            tauA(a) = NaN;
        end
    else
        tauA(a) = NaN;
    end
end
yl = get(gca,'YLim');
ylim([0 yl(2)])

maxG = maxG/max(maxG);
try
    ampHalf = interp1(maxG,amp,0.5);
catch
    try
        ampHalf = interp1(maxG(maxG > 1e-5 & maxG < 1),amp(maxG > 1e-5 & maxG < 1),0.5); % catch if constant values at beginning or end
    catch
        ampHalf = NaN;
    end
end
figure;plot(amp,maxG)
hold on;
scatter(ampHalf,0.5,'xr')
xlabel('Voltage [mV]')
ylabel('Norm. gmax')
xlim([amp(1) amp(end)])
ylim([0 1])
FontResizer
FigureResizer(height,width)
tprint(fullfile(outputFolder,sprintf('ActivCurve_%s_%s_%+.3g',mcondChannel,region,ampHalf)),'-pdf-R')

tauA(maxG <= 0.0005) = NaN;    % delete nonsense fittings
figure;plot(amp,tauA)
xlabel('Voltage [mV]')
ylabel('\tau_a_c_t [ms]')
xlim([amp(1) amp(end)])
yl = get(gca,'YLim');
ylim([0 10^ceil(log10(yl(2)))/(1+(yl(2)<10^ceil(log10(yl(2)))/2)+2*(yl(2)<10^ceil(log10(yl(2)))/4))])  % round ymax to next order of magnitude (and half/quarter it if necessary)
% title('Activation time constant')
FontResizer
FigureResizer(height,width)
tprint(fullfile(outputFolder,sprintf('ActivTime_%s_%s',mcondChannel,region)),'-pdf-R')


%% inactivation
if ~isempty(strfind(options,'-h'))
    holding_voltage = [-150,-150];
    amp = -150:5:0;
else
    holding_voltage = [-120,50];
    amp = -120:5:50;
end
if ~isempty(strfind(options,'-si'))
    fac = double(cell2mat(textscan(options,'-si%d')));
    if isempty(fac)
        fac = 1;
    end
    neuronTest.params.dt = 0.25*fac;
    predur = 50;
    dur = 10000*fac;
elseif ~isempty(strfind(options,'-fi'))
    fac = double(cell2mat(textscan(options,'-fi%d')));
    if isempty(fac)
        fac = 1;
    end
    neuronTest.params.dt = 0.05/fac;
    predur = 10;
    dur = 100/fac;
else
    neuronTest.params.dt = 0.025;
    predur = 50;
    dur = 1000;
end
if ~isempty(strfind(options,'-sa'))
    dur = dur * 5 * afac;
end
[~,out] = t2n_voltSteps(neuronTest,tree,amp,[predur dur actdur],holding_voltage);
maxG = zeros(numel(amp),1);
tauI = maxG;

figure;
hold all
xlabel('Time [ms]')
ylabel('Conductance of channel')
for a = 1:numel(amp)
    thisVec = out{a}.record{1}.cell.(cond_channel){1};                          % conductance vector of current amp
    maxG(a) = max(thisVec(out{a}.t>=predur+dur & out{a}.t<=predur+actdur+dur)); % maximal conductance during final activation
    [maxG1,indStart] = max(thisVec(out{a}.t>=predur & out{a}.t<predur+actdur));% maximal conductance during inactivation
    if indStart > 1
        indStart = indStart + find(out{a}.t>=predur,1,'first') -1;                     % correct indStart
        indEnd = find(out{a}.t>=predur+dur,1,'first')-1;                            % get end index of inactivation
        indEnd2 = find(thisVec(indStart:indEnd) < (maxG1-thisVec(indEnd))/200+thisVec(indEnd),1,'first') + indStart -1;% get value at which G during inactivation becomes smaller than 1/50 of max G
        if isempty(indEnd2)
            indEnd2 = indEnd;
        end
        lastind = find(thisVec(indStart:indEnd2)>thisVec(indEnd2),1,'last') + indStart -1; % find last point which is greater than value to be subtracted (in order to have nonzeros when applying logarithmus)
        % get part of vector which is linear, i.e. calculate derivative
        % and check where it deviates from the median by more than one
        % s.e.m.
        firstDiff = diff(log(thisVec(indStart:lastind)-thisVec(indEnd2)),1,1);
        indStartNew = find(abs(firstDiff-median(firstDiff))<=(std(firstDiff)/sqrt(numel(firstDiff))*1.01),1,'first')+indStart;   %1.01 because of float precision
        indEndNew = find(abs(firstDiff-median(firstDiff))<=(std(firstDiff)/sqrt(numel(firstDiff))*1.01),1,'last')+indStart+1;       %1.01 because of float precision
        if numel(indStartNew:indEndNew) > 1
            [ft,~] = polyfit(out{a}.t(indStartNew:indEndNew),log(thisVec(indStartNew:indEndNew)-thisVec(indEnd2)),1);
            if imag(ft(1)) ~= 0
                ft = [NaN NaN];
            end
            plot(out{a}.t,thisVec)
            plot(out{a}.t(indStartNew:indEndNew),exp(ft(1)*out{a}.t(indStartNew:indEndNew)+ft(2))+thisVec(indEnd2),'r--')
            tauI(a) = -1/ft(1);
        else
            tauI(a) = NaN;
        end
    else
        tauI(a) = NaN;
    end
end
yl = get(gca,'YLim');
ylim([0 yl(2)])

maxG = maxG/max(maxG);
if maxG(end) > 0.5 && ~any(tauI < 100)
    warning('It seems that the channel has either no or a very slow inactivation time constant. In the latter case, try to rerun this function with the option -si')
end
try
    ampHalf = interp1(maxG,amp,0.5);
catch
    ampHalf = NaN;
end
figure;plot(amp,maxG)
hold on;
scatter(ampHalf,0.5,'xr')
xlabel('Voltage [mV]')
ylabel('Norm. gmax')
xlim([amp(1) amp(end)])
ylim([0 1])
FontResizer
FigureResizer(height,width)
tprint(fullfile(outputFolder,sprintf('InactivCurve_%s_%s_%+.3g',mcondChannel,region,ampHalf)),'-pdf-R')

tauI(maxG >= 0.99) = NaN;   % delete nonsense fittings
figure;plot(amp,tauI)
xlabel('Voltage [mV]')
ylabel('\tau_i_n_a_c_t [ms]')
xlim([amp(1) amp(end)])
yl = get(gca,'YLim');
ylim([0 10^ceil(log10(yl(2)))/(1+(yl(2)<10^ceil(log10(yl(2)))/2)+2*(yl(2)<10^ceil(log10(yl(2)))/4))])  % round ymax to next order of magnitude (and half it if necessary)
% title('Inactivation time constant')
FontResizer
FigureResizer(height,width)
tprint(fullfile(outputFolder,sprintf('InactivTime_%s_%s',mcondChannel,region)),'-pdf-R')

