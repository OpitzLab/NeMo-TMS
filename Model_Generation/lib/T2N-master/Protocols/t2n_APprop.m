function [props, fig, figname] = t2n_APprop(targetfolder_data,neuron,currstep,tree,dataset,askfile)

% This function analyzes and optionally plots various action potential (AP) 
% characteristics at a current injection step previously simulated with t2n_currSteps
%
% INPUTS
% targetfolder_data destination of temporary results file
% neuron            t2n neuron structure with already defined mechanisms
% currstep          current injection amplitude at which AP properties
%                   should be analyzed [nA]
% tree              (optional) tree cell array with morphologies to
%                   calculate dependencies on tree surfaces
% dataset           (optional, so far only compatible with GC model) dataset of
%                   experimental data which should be loaded
% askfile           (optional) boolean if it should be ask for which
%                   simulation should be loaded and analyzed
%
% OUTPUTS
% props             structure with AP properties. If this is the only
%                   output, figures are not plotted
% fig               figure handles to figures
% figname           name of each figure
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

ap = -10; % minimum AP amplitude threshold for detection [mV]
if ~exist('askfile','var')
    askfile = 0;
end
newrate = 0.005; %ms
if ~exist('currstep','var')
    currstep = 0.09;
end
if ~exist('tree','var')
    tree = [];
end

if exist('dataset','var')
    uu = 1:2;
else
    uu = 2;
end
if any(uu==2)
    if askfile
        [fname,pname] = uigetfile('.mat','select Model Exp_Spiking file',t2n_catName(targetfolder_data,'Exp_Spiking',neuron.experiment,'.mat'));
        if isempty(fname) || isnumeric(fname) && fname == 0
            uu(uu==2) = [];
        else
            sim1 = load(fullfile(pname,fname));
        end
    else
        sim1 = load(t2n_catName(targetfolder_data,'Exp_Spiking',neuron.experiment,'.mat'));
    end
end

if any(uu == 1)
    [exp_iclamp,cstepsSpiking,rate] = load_ephys(dataset,'CClamp');
end
% find step which should be analyzed
s2=find(sim1.cstepsSpikingModel == currstep);
if any(uu == 1)
    s1=find(cstepsSpiking == currstep);
end


dVthresh = 15;

ylab = {'Spike width [ms]','Spike threshold [mV]','fAHP [mV]','Interspike interval [ms]','Spike amplitude [mV]','Spike width [ms]','fAHP [mV]'};%sprintf('Spike threshold [mV] (dV > %g mV/ms)',dVthresh)
figname = {'APwidth','APvthresh','APfahp','APisi','APamp','APvthreshVSAPwidth','APampVSfahp','APampabsVSfahpabs','APfahpdeact','SurfacesAPwidth','SurfacesAPnum','SurfacesfAHPabs','SurfacesCthresh','SurfacesAPamp'};

if nargout ~= 1
    for f = 1:9+~isempty(tree)*5
        fig(f) = figure;hold all,
        if f <= 5
            if f == 4
                xlabel('After n''th spike')
            else
                xlabel('n''th spike')
            end
            ylabel(ylab{f})
        else
            switch f
                case 6
                    xlabel('Spike threshold [mV]')
                    ylabel('fAHP [mV]')
                case 7
                    xlabel('Spike amplitude [mV]')
                    ylabel('Spike width [ms]')
                case 8
                    xlabel('Absolute spike amplitude [mV]')
                    ylabel('fAHP potential [mV]')
                case 9
                    xlabel('n''th spike')
                    ylabel('fAHP deactivation tau [ms]')
                case {10,12,14}
                    xlabel('surface soma [µm²]')
                case {11,13}
                    xlabel('surface dendrite [µm²]')
            end
        end
    end
end
if numel(uu) == 1
    if uu == 1
        props.APt = cell(1,size(exp_iclamp,2));
    else
        props.APt = cell(1,size(sim1.voltVec,1));
    end
else
    props.APt = cell(2,max([size(exp_iclamp,2),size(sim1.voltVec,1)]));
end
markerstyle = {'x','o','d'};
props.APiv = props.APt; props.APit = props.APt; props.APwidth = props.APt; props.APind = props.APt; props.APamp = props.APt;props.APampabs = props.APt;props.APISI = props.APt;props.fAHP = props.APt;props.fAHPabs = props.APt;props.APdeact = props.APt;
props.APic = cell(numel(uu),1);props.maxDV = props.APt;
for u = 1:numel(uu)
    switch uu(u)
        case 1
            if isempty(s1)
                continue
            else
                thisv = exp_iclamp(:,:,s1);
                thist = repmat((1/rate:1/rate:size(thisv,1)/rate)',1,size(thisv,2));
                thisv = interp1(thist(:,1),thisv,thist(1,1):newrate:thist(end,1));
                thist = repmat((thist(1,1):newrate:thist(end,1))',1,size(thisv,2));
                spikdetect = squeeze(any(exp_iclamp > -30,1));
                csteps = cstepsSpiking;
            end
        case 2
            if isempty(s2)
                continue
            end
            thisv = cell2mat(sim1.voltVec(:,s2)');
            thist = cell2mat(sim1.timeVec(:,s2)');
            thisv = interp1(thist(:,1),thisv,thist(1,1):newrate:thist(end,1));
            if size(thisv,1) == 1
                thisv = thisv';
            end
            thist = repmat((thist(1,1):newrate:thist(end,1))',1,size(thisv,2));
            spikdetect = cellfun(@(x) any(x > -30),sim1.voltVec);
            csteps = sim1.cstepsSpikingModel;
            if ~isempty(tree)
                for t = 1:numel(tree)
                    su=surf_tree(tree{t});
%                     len=len_tree(tree{t});
%                     B = B_tree(tree{t});
%                     Bs(t) = sum(B(tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')) & tree{t}.R ~= find(strcmp(tree{t}.rnames,'soma'))));
%                     lendend(t) = sum(len(tree{t}.R ~= find(strcmp(tree{t}.rnames,'soma')) & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon'))));
                    sudend(t) = sum(su(tree{t}.R ~= find(strcmp(tree{t}.rnames,'soma')) & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon'))));
                    susom(t) = sum(su(tree{t}.R == find(strcmp(tree{t}.rnames,'soma'))));
%                     Dsom(t) = tree{t}.D(1);
%                     suAIS(t) = sum(su(tree{t}.R == find(strcmp(tree{t}.rnames,'axonh'))));
                end
            end
%         case 3
%             thisv = cell2mat(sim2.voltVec(:,s2)');
%             thist = cell2mat(sim2.timeVec(:,s2)');
%             thisv = interp1(thist(:,1),thisv,thist(1,1):newrate:thist(end,1));
%             if size(thisv,1) == 1
%                 thisv = thisv';
%             end
%             thist = repmat((thist(1,1):newrate:thist(end,1))',1,size(thisv,2));
%             spikdetect = cellfun(@(x) any(x > -30),sim2.voltVec);
%             csteps = sim2.cstepsSpikingModel;
    end
    [~,props.APic{u}] = find(spikdetect == 1 & cumsum(spikdetect,2) == 1);
    props.APic{u} = csteps(props.APic{u})*1000;
    thisdv = diff(thisv,1,1)./diff(thist,1,1);
    
    for t = 1:size(thisv,2)
        flag = false;
        flag2 = false;
        [props.APampabs{u,t},props.APt{u,t}] = findpeaks(thisv(:,t),thist(:,t),'MinPeakHeight',ap,'MinPeakDistance',2); % spikes need to be 2 ms apart at least
        [~,~,props.APind{u,t}] = intersect(props.APt{u,t},thist(:,t));
        
        props.APt{u,t} = thist(props.APind{u,t});
        props.APISI{u,t} = diff(props.APt{u,t});
        count = 1;
        for n = 1 : numel(thisdv(:,t))  % go through voltage trace
            if  numel(props.APind{u,t}) >= count && thisdv(n,t) > 0 && (numel(props.maxDV{u,t}) < count || props.maxDV{u,t}(count) < thisdv(n,t))  % get time and value and voltage level of maximum voltage change
                props.maxDV{u,t}(count) = thisdv(n,t);
                props.APiv2{u,t}(count) = thisv(n,t);
                props.APit2{u,t}(count) = thist(n);
            end
            
            if numel(props.maxDV{u,t}) == count && n > props.APind{u,t}(count) && ~flag2 && thisdv(n,t)<0 && thisv(n,t) <= props.APiv2{u,t}(count)
                props.APwidth2{t}(count) = thist(n) - props.APit2{u,t}(count);  % AP width at level of fast rising phase
                flag2 = true;
            end
            
            if ~flag && thisdv(n,t) >= dVthresh && numel(props.APt{u,t})>=count && (3 > props.APt{u,t}(count)-thist(n)) %thisdv(n,t) >= dv_base(t)+2*stddv_base(t) %&&
                props.APit{u,t}(count) = thist(n);
                props.APiv{u,t}(count) = thisv(n,t);
                flag = true;
            elseif flag && thisv(n,t) <= props.APiv{u,t}(count) % spike is over
                
                props.APamp{u,t}(count) = props.APampabs{u,t}(count) - props.APiv{u,t}(count);
                props.APwidth{u,t}(count) = thist(props.APind{u,t}(count)+find(thisv(props.APind{u,t}(count)+1:end,t) <= props.APiv{u,t}(count)+props.APamp{u,t}(count)/2,1,'first'),t)-thist(find(thisv(1:props.APind{u,t}(count),t) <= props.APiv{u,t}(count)+props.APamp{u,t}(count)/2,1,'last'),t);  %width at half amplitude (as Brenner 2005)
                if props.APt{u,t}(count) < 255.5 - 5  % security zone of 5 ms not to measure props.fAHP in spikes at end of stimulation
                    if numel(props.APt{u,t}) <= count
                        tlim = 255.5;         % set window to search for props.fAHP, 255 is end of stimulation
                    else
                        tlim = min(255.5,props.APt{u,t}(count+1));
                    end
                    ind = find(thist(:,t) >= props.APt{u,t}(count) & thist(:,t) <= tlim); % get time window between spike max and next spike
%                     [props.fAHPabs{u,t}(count), ind2] = min(thisv(ind,t));  % find minimum in V in that window
%                     figure;plot(thist(:,t),thisv(:,t)),hold all;plot(thist(ind,t),thisv(ind,t))
                    [y,x] = findpeaks(-thisv(ind,t),thist(ind,t),'MinPeakDistance',1.5);
                    tim = x(find(diff(y,1)<0,1,'first'));
                    if isempty(tim)
                        [~,tim] = max(y);
                        tim = x(tim);
                    end
                    if isempty(tim)
                        props.APdeact{u,t}(count) = NaN;
                        props.fAHP{u,t}(count) = NaN;
                        props.fAHPabs{u,t}(count) = NaN;
                    else
                        props.fAHPabs{u,t}(count) = thisv(thist(:,t)==tim,t);
                        if tim - props.APt{u,t}(count) > 5 %(thist(ind(ind2),t) - props.APt{u,t}(count)) > 5 % if time from spike to fAHP is longer than 5ms then it is propably no fAHP
                            props.APdeact{u,t}(count) = NaN;
                            [~,ind2] = min(diff(thisv(ind,t),1,1));  % find minimum of derivative
                            ind2 = ind(ind2-1+find(diff(thisv(ind(ind2:end),t),1,1)> -0.005,1,'first'));  %attempt to rescue fAHP calculation..check were derivative negative steepness is over
                            if (isempty(ind2) && thist(ind(end)) == 255.5) || thist(ind2,t) - props.APt{u,t}(count) > 5
                                props.fAHP{u,t}(count) = NaN;
                                props.fAHPabs{u,t}(count) = NaN;
                            else
                                props.fAHPabs{u,t}(count) = thisv(ind2,t);
                                props.fAHP{u,t}(count) = props.APiv{u,t}(count)-props.fAHPabs{u,t}(count);
                                display('fAHP rescued')
                            end
                        else
                            props.fAHP{u,t}(count) = props.APiv{u,t}(count)-props.fAHPabs{u,t}(count);
                            
                            thist2 = thist(ind(thist(ind,t)==tim):end,t);
                            thist2 = thist2(1:find(thist2(1)+5 <= thist2,1,'first')); % nur bis +5ms
                            
                            thisv2 = thisv(ind(thist(ind,t)==tim):ind(thist(ind,t)==tim)+numel(thist2)-1,t);
                            
                            lthisv2 = log(thisv2(end)-thisv2);
                            imagind = find(imag(lthisv2)~=0 | isinf(lthisv2),1,'first')-1;
                            %                     min(imagind,find(thist2-thist(ind(ind2)-3)
                            if imagind < 4
                                props.APdeact{u,t}(count) = NaN;
                            else
                                [p,~] = polyfit(thist2(1:imagind)-thist(ind(thist(ind,t)==tim)) , lthisv2(1:imagind),1);
                                props.APdeact{u,t}(count) = -1/p(1);
                            end
                        end
                    end
                end
                flag = false;
                flag2 = false;
                count = count + 1;
            end
        end
        switch uu(u)
            case 1
                col = [0 0 0];
            case 2
                if isfield(sim1.tree{t},'col')
                    col = sim1.tree{t}.col{1};
                else
                    col = [1 0 0];
                end
        end
        
        if nargout ~= 1
            figure(fig(1)),hold on
            plot(props.APwidth{u,t},'LineWidth',1,'Color',col,'Marker',markerstyle{u})
            figure(fig(2)),hold on
            plot(props.APiv{u,t},'LineWidth',1,'Color',col,'Marker',markerstyle{u})
            figure(fig(3)),hold on
            plot(props.fAHPabs{u,t},'LineWidth',1,'Color',col,'Marker',markerstyle{u})
            figure(fig(4)),hold on
            plot(props.APISI{u,t},'LineWidth',1,'Color',col,'Marker',markerstyle{u})
            figure(fig(5)),hold on
            plot(props.APamp{u,t},'LineWidth',1,'Color',col,'Marker',markerstyle{u})
%             if u == 3 && t == size(thisv,2)
%                 fields = {'APiv','fAHP';'APamp','APwidth';'APampabs','fAHPabs'};
%                 for ff = 1:size(fields,1)
%                     
%                     tmp1 = NaN(size(thisv,2),max(cellfun(@numel,props.(fields{ff,1})(u,:))));
%                     tmp2 = NaN(size(thisv,2),max(cellfun(@numel,props.(fields{ff,2})(u,:))));
%                     for tt = 1:size(thisv,2)
%                         tmp1(tt,1:numel(props.(fields{ff,1}){u,tt})) = props.(fields{ff,1}){u,tt};
%                         tmp2(tt,1:numel(props.(fields{ff,2}){u,tt})) = props.(fields{ff,2}){u,tt};
%                     end
%                     figure(fig(5+ff)),hold on
%                     p = plot(nanmean(tmp1,1),nanmean(tmp2,1));
%                     errbar(p,cat(3,nanstd(tmp2,[],1),nanstd(tmp1,[],1)))
%                 end
%                 %                 plot(props.APiv{u,t}(1:min(numel(props.APiv{u,t}),numel(props.fAHP{u,t}))),props.fAHP{u,t}(1:min(numel(props.APiv{u,t}),numel(props.fAHP{u,t}))),'LineWidth',u,'Color',col,'Marker','x')
%                 %                 figure(fig(7)),hold on
%                 %                 plot(props.APamp{u,t}(1:min(numel(props.APamp{u,t}),numel(props.APwidth{u,t}))),props.APwidth{u,t}(1:min(numel(props.APamp{u,t}),numel(props.APwidth{u,t}))),'LineWidth',u,'Color',col,'Marker','x')
%                 %                 figure(fig(8)),hold on
%                 %                 plot(props.APampabs{u,t}(1:min(numel(props.APampabs{u,t}),numel(props.fAHPabs{u,t}))),props.fAHPabs{u,t}(1:min(numel(props.APampabs{u,t}),numel(props.fAHPabs{u,t}))),'LineWidth',u,'Color',col,'Marker','x')
%                 figure(fig(9)),hold on
%                 plot(props.APdeact{u,t},'LineWidth',1,'Color',col,'Marker',markerstyle{u})
%             elseif u < 3
                figure(fig(6)),hold on
                plot(props.APiv{u,t}(1:min(numel(props.APiv{u,t}),numel(props.fAHP{u,t}))),props.fAHP{u,t}(1:min(numel(props.APiv{u,t}),numel(props.fAHP{u,t}))),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                figure(fig(7)),hold on
                plot(props.APamp{u,t}(1:min(numel(props.APamp{u,t}),numel(props.APwidth{u,t}))),props.APwidth{u,t}(1:min(numel(props.APamp{u,t}),numel(props.APwidth{u,t}))),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                figure(fig(8)),hold on
                plot(props.APampabs{u,t}(1:min(numel(props.APampabs{u,t}),numel(props.fAHPabs{u,t}))),props.fAHPabs{u,t}(1:min(numel(props.APampabs{u,t}),numel(props.fAHPabs{u,t}))),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                figure(fig(9)),hold on
                plot(props.APdeact{u,t},'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                if uu(u) > 1 && ~isempty(tree)
                    figure(fig(10)),hold on
                    plot(susom(t),mean(props.APwidth{u,t}),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                    ylabel('AP width [ms]')
                    xlim([0 500])
                    ylim([0 1.5])
                    figure(fig(12)),hold on
                    plot(susom(t),mean(props.fAHPabs{u,t}),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                    ylabel('fAHP [mV]')
                    xlim([0 500])
                    ylim([-70 -50])
                    
                    figure(fig(11)),hold on
                    plot(sudend(t),numel(props.APt{u,t}),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                    ylabel('Number of APs')
                    ylim([0 6])
                    xlim([4000 10000])

                    figure(fig(13)),hold on
                    plot(sudend(t),props.APic{u}(t),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                    ylabel('AP thresh [pA]')
                    xlim([4000 10000])
                    ylim([0 100])
                    figure(fig(14)),hold on
                    plot3(susom(t),sudend(t),mean(props.APamp{u,t}),'LineWidth',1,'Color',col,'Marker',markerstyle{u})
                    zlabel('AP amplitude [mV]')
                    xlim([0 500])
                    ylim([4000 10000])
                    zlim([0 200])
                end
%             end
        end
    end
    
end
   
for u = 1:numel(uu)
    if uu(u) == 1 && isempty(s1)
        continue
    end
    fprintf('AP current threshold for group %d is: %g +- %g pA (s.e.m.)\n',uu(u),mean(props.APic{u}),std(props.APic{u})/sqrt(numel(props.APic{u})) )
    fprintf('AP width for group %d is: %g +- %g ms (s.e.m.)\n',uu(u),mean(cellfun(@(x) x(1),props.APwidth(u,~cellfun(@isempty,props.APwidth(u,:))))),std(cellfun(@(x) x(1),props.APwidth(u,~cellfun(@isempty,props.APwidth(u,:)))))/sqrt(numel(props.APic{u})) )
    fprintf('AP voltage threshold for group %d is: %g +- %g mV (s.e.m.)\n',uu(u),mean(cellfun(@(x) x(1),props.APiv(u,~cellfun(@isempty,props.APiv(u,:))))),std(cellfun(@(x) x(1),props.APiv(u,~cellfun(@isempty,props.APiv(u,:)))))/sqrt(numel(props.APic{u})) )
    fprintf('AP fAHP for group %d is: %g +- %g mV (s.e.m.)\n',uu(u),mean(cellfun(@(x) x(1),props.fAHP(u,~cellfun(@isempty,props.fAHP(u,:))))),std(cellfun(@(x) x(1),props.fAHP(u,~cellfun(@isempty,props.fAHP(u,:)))))/sqrt(numel(props.APic{u})) )
    fprintf('AP interspike interval for group %d is: %g +- %g ms (s.e.m.)\n',uu(u),mean(cellfun(@(x) x(1),props.APISI(u,~cellfun(@isempty,props.APISI(u,:))))),std(cellfun(@(x) x(1),props.APISI(u,~cellfun(@isempty,props.APISI(u,:)))))/sqrt(numel(props.APic{u})) )
    fprintf('AP amplitude for group %d is: %g +- %g mV (s.e.m.)\n',uu(u),mean(cellfun(@(x) x(1),props.APamp(u,~cellfun(@isempty,props.APamp(u,:))))),std(cellfun(@(x) x(1),props.APamp(u,~cellfun(@isempty,props.APamp(u,:)))))/sqrt(numel(props.APic{u})) )
end
