clear
clc
global EEG EEG_sw par

fprintf('Loading settings\n');
%% 1) Load settings
% all values, parameters and variables are saved in par (global); information
% specifically concerning EEG data are saved in EEG (global)
EEG = load('C:\Users\Norman Sinnigen\Documents\MATLAB\DATA_TO_NORMAN\OnlineArtifactCorrection_code\chanlocs_from_Analyzer2_128Ch'); % load channel location information
par.desired_ch = transpose(1:63); % number of channels
% par.desired_ch = par.desired_ch(find(par.desired_ch~=32));
EEG.chanlocs = EEG.chanlocs(par.desired_ch);  % equals the length of channel locations for channels 1-63 in variable EEG.chanlocs
EEG.nbchan = length(par.desired_ch); % equals the length of the desired_ch (63) in variable EEG.nbchan
EEG.maxChInterp = 20; % set max number of channels to be interpolated

%--------------------------------------------------------------------------
% Variables for sliding window approach
%--------------------------------------------------------------------------
par.sw = 500; % smallWindow (should be a multiple of par.defsuperbigwindow) % 250
amount_of_data_points = 300000;
%--------------------------------------------------------------------------
% Variables for ICA
%--------------------------------------------------------------------------
% to obtain reliable ICA, incoming data samples should be at least a
% multiple k of n^2, where n is the number of channels and k may need to be
% 20 or larger. Here, given n = 63 and k = 20, data samples are around 80,000

ICA = 1; % ICA = 1 enables ICA algorithm; ICA = 0 deactivates ICA
par.numOfIC = 5; % number of independent components (IC); should range between 5-10, but not below 5
par.defsuperbigwindow = 80000;

% Variables for controlling 'whitening'
% firstEig and lastEig specify the range for eigenvalues that are retained
par.firstEig          = 1; % index of largest eigenvalues to keep
par.lastEig           = par.numOfIC; % index of the last (smallest) eigenvalue to keep

% Default values for fixed point ICA parameters
par.epsilon           = 0.001; % stopping criterion
par.maxNumIterations  = 50; % default = 1000 (maximum number of iterations)

% Loading ICA template (template of ideal eye artifact (ICA_EyeHor = horizontal eye
% artifact; ICA_EyeVert = vertical eye artifact))

% for 63 channels

par.ICA_EyeHor = [0.549930755575988;-1.23201812304577;0.715267969445481;-0.808599885425582;0.628412186036173;-0.385705350733647;0.459300740834362;-0.0718191156452517;0.317085800153237;0.0315854123363991;2.46435540980346;-2.40413535849089;1.45113423439688;-1.01875224696384;0.564093672615190;-0.153263122777944;-0.0554122939106310;0.101664275357278;0.162749083456124;0.188716541531526;0.316672742973114;-0.259734247713312;0.270753605169826;-0.0544533571202849;1.36672480378071;-1.20820766319548;0.792314082880604;-0.308702767157315;0.894218999098727;-0.454368800081265;0.187306088777470;0.266121002165053;-0.379889562546222;0.301341768122391;-0.136278609219451;0.323017662652316;0.0679053165073431;0.534624639707329;-0.848034738222811;0.726696410458979;-0.612499282035063;0.524562063085504;-0.172004019145040;0.351354923870577;0.0524142736685281;1.36647663893710;-1.57082177469961;1.07051334425020;-0.734917403010183;0.551953766252701;-0.0888565277669833;1.82732373964643;-2.31506135244093;1.97068497757544;-1.63472057804342;1.01646132252746;-0.456691942286595;0.439186404443516;-0.00695794209293895;2.79239010732697;-2.16770135622766;-0.301003900681065;0.144925456188972];
par.ICA_EyeVert = [2.98725807590548;2.88837029746993;0.718652016977733;0.463101724250311;-0.0835314234275092;-0.384049424604029;-0.347061447984219;-0.714090744221027;-0.545860211915621;-0.846330229170974;1.12446421349439;0.146956814728358;-0.149826428978673;-0.903255957704736;-0.408331883731307;-1.02486218294318;0.623559360071658;-0.201121586333930;-0.507256963911764;-0.758337213031535;0.0274498698937332;-0.116526369557765;-0.275087705880795;-0.436263660831126;0.249058418745284;-0.254399515765336;-0.242770981862192;-0.752796067459360;-0.203870745552025;-1.88015453949280;-0.771789199802491;0.547520433761270;0.457296771293735;-0.194923175031498;-0.354861591035607;-0.377414696386545;-0.567121509127773;1.66938155792299;1.65560446645658;0.0885887561108664;-0.146655410701016;-0.235018897309136;-0.557726104888613;-0.534554756325647;-0.800449934751895;1.09091113932930;0.513527187056816;-0.0790027105412391;-0.564579573366736;-0.338582057848763;-0.699141399437207;2.50925208851822;1.78708357914369;0.240111690868962;-0.589995971894522;-0.320201628989981;-1.13781520326495;-0.456041947037395;-0.985493406122689;-0.255142494301014;-1.17144979944261;2.97218803246976;-0.532095112959968];
% par.ICA_EyeHor(32,:) = [];
% par.ICA_EyeVert(32,:) = [];

% for 21 channels
% par.ICA_EyeHor = [0.549930755575988;-1.23201812304577;0.715267969445481;-0.808599885425582;0.628412186036173;-0.385705350733647;0.459300740834362;-0.0718191156452517;0.317085800153237;0.0315854123363991;2.46435540980346;-2.40413535849089;1.45113423439688;-1.01875224696384;0.564093672615190;-0.153263122777944;-0.0554122939106310;0.101664275357278;0.162749083456124;0.188716541531526;0.144925456188972];
% par.ICA_EyeVert = [2.98725807590548;2.88837029746993;0.718652016977733;0.463101724250311;-0.0835314234275092;-0.384049424604029;-0.347061447984219;-0.714090744221027;-0.545860211915621;-0.846330229170974;1.12446421349439;0.146956814728358;-0.149826428978673;-0.903255957704736;-0.408331883731307;-1.02486218294318;0.623559360071658;-0.201121586333930;-0.507256963911764;-0.758337213031535;-0.532095112959968];

%--------------------------------------------------------------------------
% Threshold for calculation of spatial correlation
%--------------------------------------------------------------------------
par.th_hor = 0.7; % threshold for c-value of horizontal eye movement
par.th_ver = 0.7; % threshold for c-value of vertical eye movement

% Initialization
var_ch_sbw = zeros(EEG.nbchan, par.defsuperbigwindow/par.sw);
var_w_sbw = zeros(EEG.nbchan, par.defsuperbigwindow/par.sw);

NaN_CORRECT_window = NaN(length(par.desired_ch),par.sw);
receive = [];
ica_time_all = [];
interpolation_all = [];

%% Start labstreaming layer and open inlet for EEG data
fprintf('Open inlet\n');
lib = lsl_loadlib();
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG'); % resolve EEG stream
end
inlet = lsl_inlet(result{1}); % opens inlet for EEG data

%% Collecting data points for Independent Component Analysis
fprintf('Collecting data points for ICA\n');

superbigwindow = zeros(EEG.nbchan, par.defsuperbigwindow);
for begin = 1 : par.defsuperbigwindow
    [EEGData,~] = inlet.pull_sample(); % ts = time stamp
    EEGData = transpose(EEGData);
%     EEGData(32,:) = [];
    superbigwindow(:,begin) = EEGData;
end
raw = superbigwindow;
NaN_CORRECT = superbigwindow;

% -------------------------------------------------------------------------
% detect bad channels
% -------------------------------------------------------------------------
endd = par.sw;
col = 1;
for begin = 1:endd:par.defsuperbigwindow
    var_ch_sbw(:,col) = var(superbigwindow(:,begin:endd),0,2);
    endd = endd + par.sw;
    col = col + 1;
end
mean_ch_sbw = mean(var_ch_sbw,2) + 1.5*std(var_ch_sbw,0,2);

% -------------------------------------------------------------------------
% correct data points of superbigwindow based on calculated threshold
% -------------------------------------------------------------------------
fprintf('Correcting superbigwindow\n');
begin = 1;
endd = par.sw;

begin_NaN = 1;
endd_NaN = par.sw;

for col = 1:par.defsuperbigwindow/par.sw
    flag0 = 0;
    flag1 = 1;
    
    % find all bad channels of each window
    bad_channel = find(var_ch_sbw(:,col) > mean_ch_sbw);
    
    %%  Interpolation
    
    % if there are no bad channels, don't interpolate
    if ~isempty(bad_channel)
        
        % only interpolate bad channels if not exceeding limit
        if length(bad_channel) > EEG.maxChInterp
            superbigwindow(:,begin:endd) = [];
            flag0 = 1;
            
            NaN_CORRECT(:,begin_NaN:endd_NaN) = NaN;
            flag1=2;
        else
            
            EEG_sw = EEG; % create EEG_sw field, based on EEG-field data for following interpolation step
            EEG_sw.data = superbigwindow(:,begin:endd);
            
            [tmpdata] = Interpolation(bad_channel, EEG, EEG_sw);
            
            superbigwindow(:,begin:endd) = tmpdata;
            NaN_CORRECT(:,begin_NaN:endd_NaN) = tmpdata;
            
        end
    end
    
    if flag1==2
        begin_NaN = begin_NaN + par.sw;
        endd_NaN = endd_NaN + par.sw;
    end
    
    if flag0 ~= 1
        endd = endd + par.sw;
        begin = begin + par.sw;
        
        begin_NaN = begin_NaN + par.sw;
        endd_NaN = endd_NaN + par.sw;
    end
    
end

% -------------------------------------------------------------------------
% append missing data points
% -------------------------------------------------------------------------
fprintf('Appending missing data to superbigwindow\n');
while size(superbigwindow,2) < par.defsuperbigwindow
    
    % get missing smallWindows
    append_data = zeros(EEG.nbchan,par.sw);
    for begin = 1 : par.sw
        [EEGData,~] = inlet.pull_sample();
        EEGData = transpose(EEGData);
%         EEGData(32,:) = [];
        append_data(:,begin) = EEGData;
    end
    raw = [raw, append_data];
    
    var_append_data = var(append_data,0,2);
    
    % find all bad channels
    bad_channel = find(var_append_data > mean_ch_sbw);
    
    %  Interpolation
    if ~isempty(bad_channel)
        if length(bad_channel) > EEG.maxChInterp
            append_data = [];
            NaN_CORRECT = [NaN_CORRECT, NaN_CORRECT_window];
        else
            EEG_sw = EEG; % create EEG_sw field, based on EEG-field data for following interpolation step
            EEG_sw.data = append_data;
            
            [tmpdata] = Interpolation(bad_channel, EEG, EEG_sw);
            
            superbigwindow = [superbigwindow, tmpdata];
            NaN_CORRECT = [NaN_CORRECT, tmpdata];
            
            endd = endd + par.sw;
            col = col + 1;
            begin = begin + par.sw;
        end
    else
        % concatenate good smallWindow to superbigwindow
        superbigwindow = [superbigwindow, append_data];
        NaN_CORRECT = [NaN_CORRECT, append_data];
        
        endd = endd + par.sw;
        col = col + 1;
        begin = begin + par.sw;
    end
end

% -------------------------------------------------------------------------
% Update thresholds
% -------------------------------------------------------------------------
fprintf('Update thresholds\n');
% bad channel threshold
var_ch_sbw = zeros(EEG.nbchan, par.defsuperbigwindow/par.sw);

endd = par.sw;
col = 1;
for begin = 1 : par.sw : par.defsuperbigwindow
    var_ch_sbw(:,col) = var(superbigwindow(:,begin:endd),0,2);
    endd = endd + par.sw;
    col = col + 1;
end
mean_ch_sbw = mean(var_ch_sbw,2) + 1.5*std(var_ch_sbw,0,2);

% bad window threshold
endd = par.sw;
col = 1;
for begin = 1 : par.sw : par.defsuperbigwindow
    var_w_sbw(:,col) = mean(var(superbigwindow(:,begin:endd), 0, 1));
    endd = endd + par.sw;
    col = col + 1;
end
mean_w_sbw = mean(var_w_sbw,1) + 2*std(var_w_sbw,0,1);
mean_mean_w_sbw = mean(mean_w_sbw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Online Correction (muscle + eye) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Start online correction\n');

while size(raw,2) < amount_of_data_points
receive_tic = tic;
    superbigwindow = superbigwindow(:,end-par.defsuperbigwindow+par.sw+1:end);
    % -----------------------------------------------------------------
    % get next small window
    % -----------------------------------------------------------------
    smallWindow = zeros(EEG.nbchan,par.sw);
    for begin_sw = 1 : par.sw
        [EEGData,~] = inlet.pull_sample(); % ts = time stamp
        EEGData = transpose(EEGData);
%         EEGData(32,:) = [];
        smallWindow(:,begin_sw) = EEGData;
    end
    raw = [raw, smallWindow];
    
    mean_var_w_sw = mean(var(smallWindow,0,1));
   
    receive_toc = toc(receive_tic);
    receive = [receive, receive_toc];

    while mean_var_w_sw > mean_mean_w_sbw
        
        NaN_CORRECT = [NaN_CORRECT, NaN_CORRECT_window];
        
        % -----------------------------------------------------------------
        % get next small window
        % -----------------------------------------------------------------
        smallWindow = zeros(EEG.nbchan,par.sw);
        for begin_sw = 1 : par.sw
            [EEGData,~] = inlet.pull_sample(); % ts = time stamp
            EEGData = transpose(EEGData);
%             EEGData(32,:) = [];
            smallWindow(:,begin_sw) = EEGData;
        end
        raw = [raw, smallWindow];
        
        mean_var_w_sw = mean(var(smallWindow,0,1));
    end
    
    superbigwindow = [superbigwindow,smallWindow];
    
    % -------------------------------------------------------------------------
    % detect bad channels
    % -------------------------------------------------------------------------
  
    var_bad_ch_sw = var(superbigwindow(:,end-par.sw+1:end),0,2);
    bad_channel = find(var_bad_ch_sw > mean_ch_sbw);
   
    % if there are no bad channels, don't interpolate
    if ~isempty(bad_channel)
        
        % only interpolate bad channels if not exceeding limit
        if length(bad_channel) > EEG.maxChInterp
            superbigwindow(:,end-par.sw+1:end) = [];
            NaN_CORRECT = [NaN_CORRECT, NaN_CORRECT_window];
            %             bad_window = 1;
            %             bad_window = [bad_window, bad_window];
            continue
        else
            
            fprintf('Interpolation\n');
            EEG_sw = EEG; % create EEG_sw field, based on EEG-field data for following interpolation step
            EEG_sw.data = smallWindow;
            
            interpolation = tic;
            [tmpdata] = Interpolation(bad_channel, EEG, EEG_sw);
            end_interpolation = toc(interpolation);
            interpolation_all = [interpolation_all, end_interpolation];
            superbigwindow(:,end-par.sw+1:end) = tmpdata;
        end
    end
    
    if ICA == 1
        fprintf('Independent Component Analysis\n');
        ICA_tic = tic;
        [superbigwindow] = ICA_calc(EEG, par, superbigwindow);
        ica_toc = toc(ICA_tic);
        ica_time_all = [ica_time_all, ica_toc];
        % Output eye-artefact-corrected smallWindow
        NaN_CORRECT = [NaN_CORRECT, superbigwindow(:,end-par.sw+1:end)];
        
        
    end
    %     bad_window = 0;
    %     bad_window = [bad_window, bad_window];
    
   
end

function [tmpdata] = Interpolation(bad_channel, EEG, EEG_sw)

EEG_sw.badchans = bad_channel;
EEG.badchans = bad_channel;

EEG_sw.goodchans = setdiff(1:EEG_sw.nbchan, EEG.badchans);
EEG_sw.oldelocs = EEG_sw.chanlocs;

EEG_sw.chanlist = { EEG.chanlocs.labels };

EEG_sw.noChannelAsCell = {};
for nochanId = 1:length(EEG_sw.badchans)
    EEG_sw.noChannelAsCell{nochanId} = EEG_sw.chanlocs(EEG_sw.badchans(nochanId)).labels;
end

EEG_sw.badchans = EEG_sw.noChannelAsCell;
EEG_sw.channel = EEG_sw.chanlist;
EEG_sw.channel = sort(setdiff(lower(EEG_sw.channel), lower(EEG_sw.badchans) ));

% Decoding channels
EEG_sw.chaninds = [];
EEG_sw.alllabs = lower({ EEG_sw.chanlocs.labels });
EEG_sw.channel = lower(EEG_sw.channel);
for ind = 1:length(EEG_sw.channel)
    EEG_sw.indmatch = find(strcmp(EEG_sw.alllabs,EEG_sw.channel{ind}));
    if ~isempty(EEG_sw.indmatch)
        for tmpi = 1:length(EEG_sw.indmatch)
            EEG_sw.chaninds(end+1) = EEG_sw.indmatch(tmpi);
        end
    else
    end
end

EEG_sw.chaninds = sort(EEG_sw.chaninds);
EEG_sw.channel = EEG_sw.chaninds;

fprintf('Removing %d channel(s)...\n', EEG_sw.nbchan - length(EEG_sw.channel) );

% Performing removal
EEG_sw.diff1 = setdiff(1:size(EEG_sw.data,1), EEG_sw.channel);
EEG_sw.data(EEG_sw.diff1, :) = [];

EEG_sw.pnts = size(EEG_sw.data,2);
EEG_sw.nbchan = length(EEG_sw.channel);
EEG_sw.chanlocs = EEG_sw.chanlocs(EEG_sw.channel);
EEG_sw.chanlocs = EEG_sw.oldelocs ; % re-save chanlocs

fprintf('Interpolating missing (%d) channel(s)...\n', length(EEG.badchans));

% Find non-empty good channels and update some variables
EEG_sw.origoodchans = EEG_sw.goodchans;
EEG_sw.chanlocs = EEG.chanlocs;
EEG_sw.nonemptychans = find(~cellfun('isempty', {EEG_sw.chanlocs.theta}));

[EEG_sw.tmp, EEG_sw.indgood] = intersect(EEG_sw.goodchans, EEG_sw.nonemptychans);
EEG_sw.indgood = EEG_sw.indgood';
EEG_sw.goodchans = EEG_sw.goodchans(sort(EEG_sw.indgood));

% Getting data channel
EEG_sw.datachans = EEG_sw.goodchans;
EEG.badchans = sort(EEG.badchans);

for index_interp = length(EEG.badchans):-1:1
    EEG_sw.datachans(EEG_sw.datachans > EEG.badchans(index_interp)) = EEG_sw.datachans(EEG_sw.datachans > EEG.badchans(index_interp))-1;
end
EEG.badchans = intersect(EEG.badchans, EEG_sw.nonemptychans);

% Scan data points: get theta, rad of electrodes
EEG_sw.tmpgoodlocs     = EEG_sw.chanlocs(EEG_sw.goodchans);
EEG_sw.xelec           = [ EEG_sw.tmpgoodlocs.X ];
EEG_sw.yelec           = [ EEG_sw.tmpgoodlocs.Y ];
EEG_sw.zelec           = [ EEG_sw.tmpgoodlocs.Z ];
EEG_sw.rad             = sqrt(EEG_sw.xelec.^2+EEG_sw.yelec.^2+EEG_sw.zelec.^2);
EEG_sw.xelec           = EEG_sw.xelec./EEG_sw.rad;
EEG_sw.yelec           = EEG_sw.yelec./EEG_sw.rad;
EEG_sw.zelec           = EEG_sw.zelec./EEG_sw.rad;
EEG_sw.tmpbadlocs      = EEG_sw.chanlocs(EEG.badchans);
EEG_sw.xbad            = [ EEG_sw.tmpbadlocs.X ];
EEG_sw.ybad            = [ EEG_sw.tmpbadlocs.Y ];
EEG_sw.zbad            = [ EEG_sw.tmpbadlocs.Z ];
EEG_sw.rad             = sqrt(EEG_sw.xbad.^2+EEG_sw.ybad.^2+EEG_sw.zbad .^2);
EEG_sw.xbad            = EEG_sw.xbad./EEG_sw.rad;
EEG_sw.ybad            = EEG_sw.ybad./EEG_sw.rad;
EEG_sw.zbad            = EEG_sw.zbad ./EEG_sw.rad;

% Spherical spline interpolation
EEG_sw.values = EEG_sw.data(EEG_sw.datachans,:);
EEG_sw.newchans = length(EEG_sw.xbad);
EEG_sw.numpoints = size(EEG_sw.values,2);

% Compute g function
% Gelec
EEG_sw.unitmat = ones(length(EEG_sw.xelec(:)),length(EEG_sw.xelec));
EEG_sw.EI = EEG_sw.unitmat - sqrt((repmat(EEG_sw.xelec(:),1,length(EEG_sw.xelec)) - repmat(EEG_sw.xelec,length(EEG_sw.xelec(:)),1)).^2 +...
    (repmat(EEG_sw.yelec(:),1,length(EEG_sw.xelec)) - repmat(EEG_sw.yelec,length(EEG_sw.xelec(:)),1)).^2 +...
    (repmat(EEG_sw.zelec(:),1,length(EEG_sw.xelec)) - repmat(EEG_sw.zelec,length(EEG_sw.xelec(:)),1)).^2);

EEG_sw.g = zeros(length(EEG_sw.xelec(:)),length(EEG_sw.xelec));
EEG_sw.m = 4; % 3 is linear, 4 is best according to Perrin's curve

for n = 1:7
    EEG_sw.L = legendre(n,EEG_sw.EI);
    EEG_sw.g = EEG_sw.g + ((2*n+1)/(n^EEG_sw.m*(n+1)^EEG_sw.m))*squeeze(EEG_sw.L(1,:,:));
end
EEG_sw.Gelec = EEG_sw.g/(4*pi);

% G-function of spherical splines
EEG_sw.unitmat = ones(length(EEG_sw.xbad(:)),length(EEG_sw.xelec));
EEG_sw.EI = EEG_sw.unitmat - sqrt((repmat(EEG_sw.xbad(:),1,length(EEG_sw.xelec)) - repmat(EEG_sw.xelec,length(EEG_sw.xbad(:)),1)).^2 +...
    (repmat(EEG_sw.ybad(:),1,length(EEG_sw.xelec)) - repmat(EEG_sw.yelec,length(EEG_sw.xbad(:)),1)).^2 +...
    (repmat(EEG_sw.zbad(:),1,length(EEG_sw.xelec)) - repmat(EEG_sw.zelec,length(EEG_sw.xbad(:)),1)).^2);

EEG_sw.g = zeros(length(EEG_sw.xbad(:)),length(EEG_sw.xelec));

for n = 1:7
    EEG_sw.L = legendre(n,EEG_sw.EI);
    EEG_sw.g = EEG_sw.g + ((2*n+1)/(n^EEG_sw.m*(n+1)^EEG_sw.m))*squeeze(EEG_sw.L(1,:,:));
end
EEG_sw.Gsph = EEG_sw.g/(4*pi);

% Compute solution for parameter C
EEG_sw.meanvalues = mean(EEG_sw.values);
EEG_sw.values = EEG_sw.values - repmat(EEG_sw.meanvalues, [size(EEG_sw.values,1) 1]); % make mean zero
EEG_sw.values = [EEG_sw.values; zeros(1,EEG_sw.numpoints )];
EEG_sw.C = pinv([EEG_sw.Gelec; ones(1,length(EEG_sw.Gelec))]) * EEG_sw.values;
clear EEG_sw.values;
EEG_sw.allres = zeros(EEG_sw.newchans, EEG_sw.numpoints );

% Apply results
for j = 1:size(EEG_sw.Gsph,1)
    EEG_sw.allres(j,:) = sum(EEG_sw.C.*repmat(EEG_sw.Gsph(j,:)',[1 size(EEG_sw.C,2)]));
end

EEG_sw.allres = EEG_sw.allres + repmat(EEG_sw.meanvalues, [size(EEG_sw.allres,1) 1]);
EEG_sw.badchansdata = EEG_sw.allres;
EEG_sw.tmpdata = zeros(length(EEG.badchans), EEG_sw.pnts);
EEG_sw.tmpdata (EEG_sw.origoodchans, :) = EEG_sw.data;
EEG_sw.tmpdata (EEG.badchans,:) = EEG_sw.badchansdata;
tmpdata = EEG_sw.tmpdata;
end

function [superbigwindow] = ICA_calc(EEG, par, superbigwindow)

% In this section ICA decomposition is performed by using fastICA
% (Hyparinen, A.; 1999). fastICA seeks an orthogonal rotation of prewhitened
% data, through a fixed-point iteration scheme, that maximizes a measure
% of non-Gaussianity of the rotated components.
% fastICA ouputs the estimated separating matrix W and the corresponding
% mixing matrix A.

% Remove current EEG ICA weights and spheres
EEG_sbw.icaweights = [];
EEG_sbw.icasphere  = [];
EEG_sbw.icawinv    = [];
EEG_sbw.icaact     = [];

% Data must be in double precision
EEG_sbw.tmpdata     = superbigwindow;
EEG_sbw.tmpdata     = EEG_sbw.tmpdata  - repmat(mean(EEG_sbw.tmpdata,2), [1 size(EEG_sbw.tmpdata,2)]); % zero mean

% Begin fastICA by removing the

[EEG_sbw.tmpdata, EEG_sbw.mixedmean] = remmean(EEG_sbw.tmpdata);

%------------------------------
% Calculate PCA
%------------------------------

% Calculates the PCA matrices for given data (row) vectors. Returns
% the eigenvector (E) and diagonal eigenvalue (D) matrices containing the
% selected subspaces. Dimensionality reduction is controlled with
% the parameters 'firstEig' and 'lastEig'.

% Calculate the eigenvalues and eigenvectors of covariance
% matrix.
if fprintf ('Calculating covariance...\n'); end

[EEG_sbw.E, EEG_sbw.D] = eig(cov(EEG_sbw.tmpdata', 1));

% Sort the eigenvalues - decending.
% eigenvalues = flipud(sort(diag(EEG_sbw.D)));
eigenvalues = sort(diag(EEG_sbw.D), 'descend');


% Drop the smaller eigenvalues
if par.lastEig < EEG.nbchan
    lowerLimitValue = (eigenvalues(par.lastEig) + eigenvalues(par.lastEig + 1)) / 2;
else
    lowerLimitValue = eigenvalues(EEG.nbchan) - 1;
end

lowerColumns = diag(EEG_sbw.D) > lowerLimitValue;

% Drop the larger eigenvalues
if par.firstEig > 1
    higherLimitValue = (eigenvalues(par.firstEig - 1) + eigenvalues(par.firstEig)) / 2;
else
    higherLimitValue = eigenvalues(1) + 1;
end
higherColumns = diag(EEG_sbw.D) < higherLimitValue;

% Combine the results from above
selectedColumns = lowerColumns & higherColumns;

% Select the colums which correspond to the desired range
% of eigenvalues (eigenvalues and eigenvectors are not sorted).

% for eigenvector E:
numTaken = 0;
for i = 1 : size (selectedColumns, 1)
    if selectedColumns(i, 1) == 1
        takingMask(1, numTaken + 1) = i; %#ok<AGROW>
        numTaken = numTaken + 1;
    end
end
EEG_sbw.E = EEG_sbw.E(:, takingMask); % --> E: Eigenvector matrix

% for diagonal eigenvalue D:
clear numTaken takingMask
numTaken = 0;
selectedColumns = selectedColumns';
EEG_sbw.D = EEG_sbw.D';
for i = 1 : size (selectedColumns, 2)
    if selectedColumns(1, i) == 1
        takingMask(1, numTaken + 1) = i;
        numTaken = numTaken + 1;
    end
end
EEG_sbw.D = EEG_sbw.D(:, takingMask)';

clear takingMask numTaken
numTaken = 0;
selectedColumns = selectedColumns';
for i = 1 : size (selectedColumns, 1)
    if selectedColumns(i, 1) == 1
        takingMask(1, numTaken + 1) = i;
        numTaken = numTaken + 1;
    end
end
EEG_sbw.D = EEG_sbw.D(:, takingMask); % --> D: Diagonal eigenvalue matrix

%-------------------
% Whitening the data
%-------------------
% The following section whitens the data (row vectors) and reduces dimension.
% Returns the whitened vectors (row vectors), whitening and dewhitening matrices.

% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
EEG_sbw.whiteningMatrix = inv(sqrt (EEG_sbw.D)) * EEG_sbw.E';
% EEG_sbw.whiteningMatrix = (sqrt (EEG_sbw.D)) / EEG_sbw.E';
EEG_sbw.dewhiteningMatrix = EEG_sbw.E * sqrt (EEG_sbw.D);

% Project to the eigenvectors of the copariance matrix.
% Whiten the samples and reduce dimension simultaneously.
% --> whitening
EEG_sbw.whitesig =  EEG_sbw.whiteningMatrix * EEG_sbw.tmpdata;

%----------------
% Calculate  ICA
%----------------

% Calculate the ICA with fixed-point algorithm
[EEG_sbw.vectorSize, EEG_sbw.numSamples]    = size(EEG_sbw.whitesig);
if fprintf('Starting ICA calculation...\n'); end

% Estimate all the independent components in parallel
EEG_sbw.A = zeros(EEG_sbw.vectorSize, par.numOfIC);  % Dewhitened basis vectors

% Take random orthonormal initial vectors
EEG_sbw.B = orth (randn (EEG_sbw.vectorSize, par.numOfIC));

EEG_sbw.BOld = zeros(size(EEG_sbw.B));
EEG_sbw.BOld2 = zeros(size(EEG_sbw.B));

% This is the actual fixed-point iteration loop

for round = 1:par.maxNumIterations + 1
    
    if round == par.maxNumIterations + 1
        fprintf('No convergence after %d steps\n', par.maxNumIterations);
        if ~isempty(EEG_sbw.B)
            % Symmetric orthogonalization
            EEG_sbw.B = EEG_sbw.B * real(inv(EEG_sbw.B' * EEG_sbw.B)^(1/2));
            
            EEG_sbw.W = EEG_sbw.B' * EEG_sbw.whiteningMatrix;
            EEG_sbw.A = EEG_sbw.dewhiteningMatrix * EEG_sbw.B;
        else
            EEG_sbw.W = [];
            EEG_sbw.A = [];
        end
        continue;
    end
    
    % Symmetric orthogonalization
    EEG_sbw.B = EEG_sbw.B * real(inv(EEG_sbw.B' * EEG_sbw.B)^(1/2));
    
    % Test for termination condition. Note that we consider opposite
    % directions here as well.
    minAbsCos = min(abs(diag(EEG_sbw.B' * EEG_sbw.BOld)));
    
    if (1 - minAbsCos < par.epsilon)
        if fprintf('Convergence after %d steps\n', round); end
        
        % Calculate the de-whitened vectors
        EEG_sbw.A = EEG_sbw.dewhiteningMatrix * EEG_sbw.B;
        break;
    end
    
    EEG_sbw.BOld2 = EEG_sbw.BOld;
    EEG_sbw.BOld = EEG_sbw.B;
    
    % Show the progress...
    if round == 1
        fprintf('Step no. %d\n', round);
    else
        fprintf('Step no. %d, change in value of estimate: %.3g \n', round, 1 - minAbsCos);
    end
    
    EEG_sbw.B = (EEG_sbw.whitesig * (( EEG_sbw.whitesig' * EEG_sbw.B) .^ 3)) / EEG_sbw.numSamples - 3 * EEG_sbw.B;
    
end

% Calculate ICA filters

EEG_sbw.W = EEG_sbw.B' * EEG_sbw.whiteningMatrix;

% ICA Output
EEG_sbw.icawinv    = EEG_sbw.A; % mixing matrix A
EEG_sbw.icaweights = EEG_sbw.W; % inverse matrix of A

% Update weight and inverse matrices etc...
if isempty(EEG_sbw.icaweights)
    EEG_sbw.icaweights = pinv(EEG_sbw.icawinv);
end
if isempty(EEG_sbw.icasphere)
    EEG_sbw.icasphere  = eye(size(EEG_sbw.icaweights,2));
end
if isempty(EEG_sbw.icawinv)
    EEG_sbw.icawinv    = pinv(EEG_sbw.icaweights*EEG_sbw.icasphere); % a priori same result as inv
end

%% 5.2) Correlation of independent components

% In this section each of the calculated component (EEG_sbw.icawinv) is correlated with
% each template (template for horizontal eye movement 'EyeHor' AND template
% for vertical eye movement 'EyeVert'). If the c (correlation coefficient)
% or p (matrix of p-values for testing the hypothesis of no correlation
% against the alternative hypothesis of a nonzero correlation) value
% exceeds the threshold (predetermined c-value) of 0.7, then reject the
% component and reconstruct the

EEG_sbw.components = []; % initialising components

for q = 1:size(EEG_sbw.icawinv,2)
    
    % Calculate c and p-value of horizontal and vertical eye movement
    
    [par.c_hor, par.p_hor] = corr(EEG_sbw.icawinv(:,q)/std(EEG_sbw.icawinv(:,q)), par.ICA_EyeHor/std(par.ICA_EyeHor));
    [par.c_vert, par.p_vert] = corr(EEG_sbw.icawinv(:,q)/std(EEG_sbw.icawinv(:,q)), par.ICA_EyeVert/std(par.ICA_EyeVert));
    
    if (abs(par.c_hor) > par.th_hor) || (abs(par.c_vert) > par.th_ver)
        EEG_sbw.components = cat(2, EEG_sbw.components, q);
    end
end

if ~isempty(EEG_sbw.components) % if there are components, reconstruct data
    
    %-----------------------
    % Reconstruction of data
    %-----------------------
    
    % Rejecting bad component
    EEG_sbw.component_keep  = setdiff_bc(1:size(EEG_sbw.icaweights,1), EEG_sbw.components);
    EEG_sbw.compproj        = EEG_sbw.icawinv(:, EEG_sbw.component_keep) * ((EEG_sbw.icaweights(EEG_sbw.component_keep,:)*EEG_sbw.icasphere)*superbigwindow);
    
    % Data reconstruction (back-projection: forward mixing
    % process from IC's to sclap channels)
    superbigwindow(par.desired_ch',:) = EEG_sbw.compproj;
    EEG_sbw.goodinds    = setdiff_bc(1:size(EEG_sbw.icaweights,1), EEG_sbw.components);
    EEG_sbw.icawinv     = EEG_sbw.icawinv(:,EEG_sbw.goodinds);
    EEG_sbw.icaweights  = EEG_sbw.icaweights(EEG_sbw.goodinds,:);
end
end