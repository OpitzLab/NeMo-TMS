function TMS_Waveform()
%% Generate TMS pulse trains and write them in a file to be used in the
% NEURON simulation.

%  Receive TMS parameters
% 1: Monophasic, 2: Biphasic
TMS_type = menu('Choose TMS pulse type:','Monophasic','Biphasic');
STEP_type = menu('Choose step time (us) :','5','25 (Default)');


%%
prompt = {'\bfEnter the inter-pulse interval in ms:',...
    '\bf Enter the number of pulses:'};
dlgtitle = 'TMS Parameters';
dims = [1 92];
opts.Interpreter = 'tex';
opts.WindowStyle = 'normal';
definput = {'1000','10'};
while 1
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    ipi = str2double(answer{1}); % inter-pulse interval
    nump = str2double(answer{2}); % number of pulses
    if ~isnan(ipi) && ~isnan(nump) && ipi>0 && nump>=1
        break
    end
    definput = answer;
    if isnan(ipi)
        definput{1} = 'Wrong format!';
    elseif ipi<=0
        definput{1} = 'Value should be positive!';
    end
    if isnan(nump)
        definput{2} = 'Wrong format!';
    elseif nump < 1
        definput{2} = 'At least one pulse is needed!';
    end
end
%% Read TMS single pulse file

if TMS_type == 1
    load(['./original_waveforms' filesep 'TMS_mono.mat']);
else
    load(['./original_waveforms' filesep 'TMS_bi.mat']);
end

% figure
% plot(TMS_t,TMS_E,'o--')
% grid minor
%% Generate pulse train
dt= 0.005;
step = 0.005;
if length(TMS_E) > round(ipi/dt)
    error('Inter-pulse interval cannot be shorter than TMS pulse duration.');
end
delay_start = 40; % delay at the beginning before TMS delivery 40
delay_end = 40; % delay after the TMS delivery 40
ipi = round(ipi/dt)*dt;
train_length = delay_start + (nump-1)*ipi + delay_end; % in ms
train_t = (0:dt:train_length)';

pulse_extend = [TMS_E; zeros(round(ipi/dt)-length(TMS_E),1)];

train_E = zeros(delay_start/dt+1,1); % zeropad before TMS delivery
train_E = [train_E; repmat(pulse_extend,nump-1,1)]; % Append TMS pulses
train_E = [train_E; TMS_E; zeros(length(train_t)-length(train_E)-length(TMS_E),1)];

%% Downsampling if dt=0.025us

if STEP_type==2
    train_E = train_E(1:5:end);
    train_t = train_t(1:5:end) ; 
    step = 0.025;
end

%% Normalization
train_E = train_E/max(train_E);

%% Save timing information
timing = [ipi;delay_start;nump;step];

%% Save train
if ~exist('../../Results/TMS_Waveform','dir')
    mkdir('../../Results/TMS_Waveform');
end
save(['../../Results/TMS_Waveform' filesep 'TMS_E_train.txt'], 'train_E','-ascii');
save(['../../Results/TMS_Waveform' filesep 'TMS_t_train.txt'], 'train_t','-ascii');
save(['../../Results/TMS_Waveform' filesep 'TMS_timing.txt'], 'timing', '-ascii');



disp('Successfully generated the TMS waveform!');
end
