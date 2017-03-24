
% Function to run spatial 3 oddity comparison of texture images. Each image
% is physically unique, and subjects compare which is odd (from two
% original / one synth or two synth / one original). Almost identical to
% Balas (2006).

function dnn_textures_expt_12()
%% Experiment number:
expt_num = 12;

%% Set random seed:
seed = sum(100*clock);
reset(RandStream.getGlobalStream,seed);

%% Ask for some starting input:
param_ok = 0;
while param_ok == 0
    in = input('Are you in the lab? [y / n]', 's');
    % check parameters:
    if strcmp(in, 'y')
        lab = true;
        use_eyetracker = true;
        param_ok = 1;
    elseif strcmp(in, 'n')
        lab = false;
        use_eyetracker = false;
        param_ok = 1;
    else
        warning('Not a valid option. Try again.');
    end
end

if use_eyetracker ~= false
    param_ok = 0;
    while param_ok == 0
        in = input('Do you want to use the eyetracker? [y / n]', 's');
        % check parameters:
        if strcmp(in, 'y')
            use_eyetracker = true;
            param_ok = 1;
        elseif strcmp(in, 'n')
            use_eyetracker = false;
            param_ok = 1;
        else
            warning('Not a valid option. Try again.');
        end
    end
end

%% Parameters

% get parameters from files:
filepathToHere=pwd;
gen_file = fullfile(filepathToHere, sprintf('params_expt_%s.yaml', num2str(expt_num)));
params = ReadYaml(gen_file);

assert(params.experiment_num == expt_num);

%% Ask for input parameters:

subj = input('Enter subject code: [practice, S1, S2,..]','s');

disp(['Check that viewing distance is ', num2str(params.dist_monitor), ' cm!!!'])

param_ok = 0;
while param_ok == 0
    ms_stim = input('Presentation time (ms)? [200, 2000]: ');
    % check parameters:
    switch ms_stim
        case 200
            param_ok = 1;
        case 2000
            param_ok = 1;
        otherwise
            warning('Not a valid option. Try again.');
            param_ok=0;
    end
end

param_ok = 0;
while param_ok == 0
    n_reps = input('N reps (default = 10)?: ');

    if ischar(n_reps)
        n_reps = num2str(n_reps);
    end

    param_ok = 1;
end

if ~lab
    param_ok = 0;
    while param_ok == 0
        in = input('fullscreen? [y / n] ', 's');
        % check parameters:
        if strcmp(in, 'y')
            fullscreen = true;
            param_ok = 1;
        elseif strcmp(in, 'n')
            fullscreen = false;
            param_ok = 1;
        else
            warning('Not a valid option. Try again.');
        end
    end
end

%% Set practice params if subject == "practice"

if strcmp(subj, 'practice')
    % practice on a few images with longer times.
    params.ms_stim = 1000;
    practice = 1;
    fixation_check = 0;
else
    practice = 0;
    fixation_check = 0;
end

%% Start the hardware
if lab
    % load gamma correction file:
    calib = load('/home/localadmin/calibration/cin_spatial_gray/cin_spatial_gray2016_06_02_1404.mat');

    clut = spline(calib.measurement, calib.input, linspace(0,1,(2^8))');  % an 8-bit luminance lut.
    %     clut(1:10) = 0;
    clut(clut<0)=0;

    % LCD-initialization
    win = window('cin_spatial_gray', 'bg_color', params.bg_color, 'clut', clut);
    listener = listener_buttonbox('names', {'Green', 'White', 'Red'});
    waiter = listener_buttonbox('names', {'Green', 'White', 'Red'}, ...
        'does_interrupt', true);
    aud_volume = 0.5;
    aud = dpixx_audio_port('volume', aud_volume);
    aud.create_beep('short_low', 'low', .15, 0.25);
    aud.create_beep('short_high', 'high', .15, 0.25);

    HideCursor(win.h);
else
    if fullscreen
        win = window('debug', 'bg_color', params.bg_color, 'rect', [0, 0, 1200, 900]);
    else
        win = window('debug', 'bg_color', params.bg_color, 'rect', [0, 0, 800, 600]);
    end
    listener = listener_keyboard('names', {'LeftArrow', 'DownArrow', 'RightArrow', 'Escape'});
    waiter = listener_keyboard('names', {'LeftArrow', 'DownArrow', 'RightArrow', 'Escape'}, ...
        'does_interrupt', true);
end


if use_eyetracker
    % initialise eyetracker
    try
        % Provide Eyelink with details about the graphics environment
        % and perform some initializations. The information is returned
        % in a structure that also contains useful defaults
        % and control codes (e.g. tracker state bit and Eyelink key values).
        el=EyelinkInitDefaults(win.h);

        % Initialization of the connection with the Eyelink Gazetracker.
        Eyelink('Initialize','PsychEyelinkDispatchCallback');

        Eyelink('Initialize','PsychEyelinkDispatchCallback');

        [~, vs]=Eyelink('GetTrackerVersion');
        fprintf('Running experiment on a ''%s'' tracker.\n', vs );

        % make sure that we get gaze data from the Eyelink
        Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');

        eye_used = -1;

    catch eyelink_error
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        % Close window:
        sca;
        % Restore keyboard output to Matlab:
        ListenChar(0);
        commandwindow;
        disp(eyelink_error);
        disp(eyelink_error.message);
        disp(eyelink_error.stack.line);
    end
end

%% Set file paths and open file to write to

data_path = fullfile(filepathToHere, sprintf('../../raw_data/experiment_%s/', num2str(params.experiment_num)));
orig_im_path = fullfile(filepathToHere, sprintf('../../stimuli/stimulus_set_4/preprocessed_ims/'));
ps_path = fullfile(filepathToHere, sprintf('../../stimuli/stimulus_set_4/ps_synths/'));
conv5_path = fullfile(filepathToHere, sprintf('../../stimuli/stimulus_set_4/conv5/'));
powerspec_path = fullfile(filepathToHere, sprintf('../../stimuli/stimulus_set_4/powerspec/'));
eye_data_path = fullfile(data_path, '/eye_data_files/');

if ~exist(data_path, 'dir')
    mkdir(data_path);
end

if ~exist(eye_data_path , 'dir')
    mkdir(eye_data_path);
end

%% Setup file saving structure
session = 1;

datafilename = fullfile(data_path,...
    sprintf('experiment-%s_sub_%s_session_%s.csv', num2str(params.experiment_num), subj, num2str(session)));

% Check for existing result file to prevent accidentally overwriting files
% from a previous session:
while fopen(datafilename, 'rt')~=-1
    fclose('all');
    warning('File already exists. Using next session number.');
    session = session +1;
    datafilename = fullfile(data_path,...
        sprintf('experiment-%s_sub_%s_session_%s.csv', num2str(params.experiment_num), subj, num2str(session)));
end

eyedata_fname = fullfile(eye_data_path, ...
    sprintf('experiment-%s_sub_%s_session_%s.edf', num2str(params.experiment_num), subj, num2str(session)));
% the remote file is stored on MS-DOS, so has filename length restrictions:
if practice
    edf_file_remote = sprintf('prac_%s.edf', num2str(session));
else
    edf_file_remote = sprintf('%s_%s.edf', subj, num2str(session));
end


%% Setup trial structure and counterbalancing

im_data_file = readtable([orig_im_path, 'stim_info.csv'], 'delimiter', ',');

im_code = table2array(im_data_file(:, 'im_code'));

% create balanced factors with n_reps:
[im_code, image_model] = BalanceFactors(n_reps, 1, im_code, {'powerspec', 'conv5', 'PS'});

% half the trials are "oddball natural":
oddball = Shuffle(repmat([0; 1], length(im_code)/2, 1));

% target position balanced:
target_vec = repmat([1; 2; 3], floor(length(im_code) / 3), 1);
while length(target_vec) < floor(length(im_code))
    target_vec = [target_vec; randi(3, 1, 1)];
end
target_loc = Shuffle(target_vec);
target_loc = num2str(target_loc);

% synth number balanced:
synth_num = repmat([1:10]', floor(length(im_code) / 10), 1);
while length(synth_num) < floor(length(im_code))
    synth_num = [synth_num; randi(10, 1, 1)];
end
synth_num = Shuffle(synth_num);

% random crop position x,y (cols, rows). half move along x, half along y.
% specifies the top left corner of the crop patch.
crop_x = [ones(floor(length(im_code) / 2), 1); ...
    randi(params.image_size_px, floor(length(im_code) / 2), 1)];
crop_y = [randi(params.image_size_px, floor(length(im_code) / 2), 1); ...
    ones(floor(length(im_code) / 2), 1)];

% random crop position two (used for second crop of same image):
crop_x_2 = [ones(floor(length(im_code) / 2), 1) * (params.image_size_px + 1); ...
    randi(params.image_size_px, floor(length(im_code) / 2), 1)];
crop_y_2 = [randi(params.image_size_px, floor(length(im_code) / 2), 1); ...
    ones(floor(length(im_code) / 2), 1) * (params.image_size_px + 1)];

% ms_stim constant for a session:
ms_stim = repmat(ms_stim, length(im_code), 1);

data = table(image_model, im_code, ...
    oddball, target_loc, ms_stim, synth_num, ...
    crop_x, crop_y, crop_x_2, crop_y_2);

n_trials = height(data); % number of trials

% fill new variables in data structure:
for trial=1:n_trials
    dat_struct(trial).subj = subj;
    dat_struct(trial).session = session;
    dat_struct(trial).test_location = 'cin';
    dat_struct(trial).rand_seed = seed;
    dat_struct(trial).trial = nan;
    dat_struct(trial).response = [];
    dat_struct(trial).rt = nan;
    dat_struct(trial).eye_valid = nan;
    dat_struct(trial).prop_eye_valid = nan;
end

% add to table:
data = [data, struct2table(dat_struct)];

%% pseudo randomization-> ensure that an image never follows itself.

% for pilot testing:
if practice
    if strcmp(subj, 'practice')
        % if practice:
        data = data(1:30,:);  % first 30 trials.
    else
        params.break_trials = 10;
    end
    n_trials = height(data); % number of trials
    trial_row_idx = Shuffle(1:n_trials); % Create random trial order
else
    trial_row_idx = pseudo_rand(data, 2); % random order; make sure images not together
end


%% Calculate spatial and temporal display workings.

% debug setup sometimes can't find framerate:
if lab == false && win.framerate == 0
    win.framerate = 60;
end

n_frames_stim = ceil(ms_stim * win.framerate / 1000);
n_frames_resp = ceil(params.ms_resp * win.framerate / 1000);
n_frames_fixate = ceil(params.ms_fixation * win.framerate / 1000);
n_frames_feedback = ceil(params.ms_feedback * win.framerate / 1000);

hann.size = [params.image_size_px, params.image_size_px];
hann.alpha = params.hann_alpha;
spatial_window = distrib_tapered_cosine(hann);

% too hard to specify these in the yaml file because cell arrays:
fix_colour_oval = [0 0 0];
fix_colour_normal = [0.7 0.7 0.7];
fix_colour_correct = [0.9 0.9 0.9];
fix_colour_wrong = [0.3 0.3 0.3];

fix_position = [win.rect(3) / 2, win.rect(4) / 2 ];

% interval angles (1 = bl, 2 = middle, 3 = br)[radians; CCW from right]:
angle_2 = 90 * (pi / 180);  % middle patch vertical
angle_1 = angle_2 + ((2 * pi) / 3);
angle_3 = angle_2 - ((2 * pi) / 3);

% load sad face texture:
fname = fullfile(filepathToHere, sprintf('../../stimuli/sadface.png'));
im = im2double(imread(fname));
sad_tex = win.make_texture(im);
sad_rect = [0, 0, params.ppd*2, params.ppd*2];
sad_rect = CenterRectOnPoint(sad_rect, fix_position(1),fix_position(2));

ok_range_px = params.ppd * params.check_dist;  % the radius of acceptable fixation, in px.

experiment_start_time = tic;

%% Calibrate eyetracker, if used.
if use_eyetracker
    try
        % open file to record data to
        Eyelink('Openfile', edf_file_remote);

        % Calibrate the eye tracker
        EyelinkDoTrackerSetup(el);

        %         % do a final check of calibration using driftcorrection
        %         EyelinkDoDriftCorrection(el);

    catch eyelink_error
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        % Close window:
        sca;
        % Restore keyboard output to Matlab:
        ListenChar(0);
        commandwindow;
        disp(eyelink_error);
        disp(eyelink_error.message);
        disp(eyelink_error.stack.line);
    end
end

%% Trials
try
    ListenChar(2);
    if IsWin()
        ListenChar(0);
    end

    for trial=1:n_trials

        this_trial_idx = trial_row_idx(trial);  % the index of this trial in the data frame.

        % Save properties in data structure
        data.trial(this_trial_idx) = trial;

        if trial == 1
            win.pause_trial(waiter, 'Press any button to start!');
            for itic = 1 : n_frames_fixate
                win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
                win.flip();
            end
        end

        %% Load parameters for this trial
        im_code = table2array(data(this_trial_idx, 'im_code'));
        image_model = table2array(data(this_trial_idx, 'image_model'));
        oddball = table2array(data(this_trial_idx, 'oddball'));
        target_loc = table2array(data(this_trial_idx, 'target_loc'));
        synth_num = table2array(data(this_trial_idx, 'synth_num'));
        crop_x = table2array(data(this_trial_idx, 'crop_x'));
        crop_y = table2array(data(this_trial_idx, 'crop_y'));
        crop_x_2 = table2array(data(this_trial_idx, 'crop_x_2'));
        crop_y_2 = table2array(data(this_trial_idx, 'crop_y_2'));

        if iscell(im_code)
            im_code = im_code{(1)};
        end

        %% target images:
        im_0_fname = get_natural_im(im_code);  % the natural image
        im_1_fname = get_model_im(im_code, image_model, synth_num);

        % determine target / nontarget texture according to oddball:
        if oddball == 0
            targ_fname = im_0_fname;
            nontarg_fname = im_1_fname;
        elseif oddball == 1
            targ_fname = im_1_fname;
            nontarg_fname = im_0_fname;
        end

        % determine interval texture according to target_pos:
        if strcmp(target_loc, '1')
            interval_1_fname = targ_fname;
            interval_2_fname = nontarg_fname;
            interval_3_fname = nontarg_fname;

            interval_1_crop = 1;
            interval_2_crop = 1;
            interval_3_crop = 2;  % second nontarg is physically diff.
        elseif strcmp(target_loc, '2')
            interval_1_fname = nontarg_fname;
            interval_2_fname = targ_fname;
            interval_3_fname = nontarg_fname;

            interval_1_crop = 1;
            interval_2_crop = 1;
            interval_3_crop = 2;  % second nontarg is physically diff.

        elseif strcmp(target_loc, '3')
            interval_1_fname = nontarg_fname;
            interval_2_fname = nontarg_fname;
            interval_3_fname = targ_fname;

            interval_1_crop = 1;
            interval_2_crop = 2;  % second nontarg is physically diff.
            interval_3_crop = 1;

        end

        % load images, make textures.
        if interval_1_crop == 1
            tex_1 = make_image_texture(interval_1_fname, crop_x, crop_y);
        elseif interval_1_crop == 2
            tex_1 = make_image_texture(interval_1_fname, crop_x_2, crop_y_2);
        end

        if interval_2_crop == 1
            tex_2 = make_image_texture(interval_2_fname, crop_x, crop_y);
        elseif interval_2_crop == 2
            tex_2 = make_image_texture(interval_2_fname, crop_x_2, crop_y_2);
        end

        if interval_3_crop == 1
            tex_3 = make_image_texture(interval_3_fname, crop_x, crop_y);
        elseif interval_3_crop == 2
            tex_3 = make_image_texture(interval_3_fname, crop_x_2, crop_y_2);
        end

        % determine drawing rects:
        stim_eccent_px = params.ppd * params.eccent_deg;  % eccent from fixation to inner edge of stimulus.
        stim_rect_base = [0, 0, params.image_size_px, params.image_size_px];

        % position of each interval:
        rect_1 = place_rect1(fix_position, stim_eccent_px,...
            angle_1, stim_rect_base, params.image_size_px);
        rect_2 = place_rect2(fix_position, stim_eccent_px, ...
            angle_2, stim_rect_base, params.image_size_px);
        rect_3 = place_rect3(fix_position, stim_eccent_px, ...
            angle_3, stim_rect_base, params.image_size_px);

        %% start eye recording, check that fixation is within 2 degrees of middle.
        % wait ITI.
        if use_eyetracker
            % start recording eye position
            Eyelink('StartRecording');
            WaitSecs(0.01);
            Eyelink('StartRecording');

            % do fixation checking if desired:
            if fixation_check
                check_fixation(params.check_dist, params.check_time, params.recal_time, ...
                    params.ppd, eye_used, fix_position, win, el)
            end

        else

        end

        if lab
            aud.play('short_high');
            win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
            win.flip();
        end

        % mark zero-plot time in data file
        if use_eyetracker
            % get ready to record eye positions:
            cum_gaze_pos = [];
            % get one sample to determine eye used:
            [cum_gaze_pos, eye_used] = cumulative_gaze_pos(cum_gaze_pos, eye_used, el);
            Eyelink('Message', ['start_trial: ' num2str(trial)]);
        end

        %% stimulus - gap intervals.

        % now we start recording button presses at the start of the
        % presentation:
        listener.start();
        % stimulus display:
        for itic = 1 : n_frames_stim
            win.draw(tex_1, 1, rect_1);
            win.draw(tex_2, 1, rect_2);
            win.draw(tex_3, 1, rect_3);
            win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
            Screen('DrawingFinished', win.h);
            if use_eyetracker
                [cum_gaze_pos, eye_used] = cumulative_gaze_pos(cum_gaze_pos, eye_used, el);
            end
            win.flip();

            %             if itic == n_frames_stim / 2
            %                 % for printing demos
            %                 imageArray = Screen('GetImage', win.h);
            %                 % save output:
            %                 out_name = fullfile(filepathToHere, ...
            %                     sprintf('expt_%s_eccent_%s_model_%s_size_%s.png', ...
            %                     num2str(params.experiment_num),...
            %                     num2str(params.eccent_deg),...
            %                     image_model{1},...
            %                     num2str(params.image_size_px)));
            %                 imwrite(imageArray, out_name);
            %             end
        end
        if use_eyetracker
            Eyelink('Message', ['end_interval: ' num2str(1)]);
        end


        % clear the screen:
        win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
        win.flip();

        % flag trial end:
        if use_eyetracker
            Eyelink('Message', 'end_trial');
        end

        if lab  % tone to signal trial end.
            aud.play('short_high');
            win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
            win.flip();
        end

        %% Wait for a response, or response interval to be exceeded.

        %%%%%% fixed timing response interval %%%%%%
        % Response interval
        for itic = 1 : n_frames_resp
            win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
            win.flip();
        end
        response = listener.stop();

        % Save the responses
        [press, rt] = listener.response.get_presses('last');

        %%%%% fixed timing response interval %%%%%%

        switch press
            case 1
                res = '1';
            case 2
                res = '2';
            case 3
                res = '3';
            otherwise
                res = 'na';
        end

        data.response(this_trial_idx) = {res};
        data.rt(this_trial_idx) = rt;

        if use_eyetracker
            Eyelink('Message', 'end_response');
            Eyelink('Message', ['response = ', res]);

            % process eye data for this trial:
            valid_vec = check_gaze_valid(cum_gaze_pos, fix_position, ok_range_px);
            [eye_valid, prop_valid] = summarise_trial_eye_valid(valid_vec);
            data.eye_valid(this_trial_idx) = eye_valid;
            data.prop_eye_valid(this_trial_idx) = prop_valid;
        end

        % close unused textures to save memory:
        Screen('Close', [tex_1, tex_2, tex_3]);

        %% Breaks and feedback

        % provide feedback via fixation cross colour ON PRACTICE TRIALS ONLY:
        if practice
            if strcmp(res, target_loc)
                for itic = 1 : n_frames_feedback
                    win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_correct);
                    win.flip();
                end
            else
                if lab
                    aud.play('short_low');
                    win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
                    win.flip();
                end
                for itic = 1 : n_frames_feedback
                    win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_wrong);
                    win.flip();
                end
            end
        end
        win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
        win.flip();

        if use_eyetracker
            Eyelink('Message', 'end_feedback');
            Eyelink('StopRecording');
            WaitSecs(0.01);
            Eyelink('StopRecording');

            % if eye was invalid, show a sad face:
            if ~eye_valid
                for itic = 1 : (n_frames_feedback * 2)
                    win.draw(sad_tex, 1, sad_rect);
                    win.flip();
                end
            end
        end

        % if this is a trial break:
        if mod(trial, params.break_trials)==0 && trial~=n_trials
            trial_sum = 0;
            for icorr=(trial-(params.break_trials-1)):trial
                this_trial_idx = trial_row_idx(icorr);
                if strcmp(data.response(this_trial_idx), data.target_loc(this_trial_idx))
                    trial_sum = trial_sum+1;
                end
            end
            perc_corr = (trial_sum/params.break_trials)*100;

            win.pause_trial(waiter, ...
                sprintf(['%d out of the last %d trials correct. '...
                'That corresponds to %.2f %% correct. \n'...
                'You have finished %d blocks out of %d. \n'...
                '\nPress a key to continue!'], ...
                trial_sum,params.break_trials,perc_corr,...
                (trial/params.break_trials), round(n_trials/params.break_trials)));

            win.draw_fixation(params.ppd, fix_position, fix_colour_oval, fix_colour_normal);
            win.flip();
            WaitSecs(params.ms_fixation/1000);
        else
            % intertrial interval
            WaitSecs(params.ms_iti/1000);
        end

        % Responded to last trial
        if trial == n_trials
            trial_sum = 0;
            for icorr=(trial-(n_trials-1)):trial
                if strcmp(data.response(icorr), data.target_loc(icorr))
                    trial_sum = trial_sum+1;
                end
            end
            perc_corr = (trial_sum/n_trials)*100;

            win.pause_trial(waiter, ...
                sprintf(['%d out of %d trials correct. '...
                'That corresponds to %.2f %% correct. \n' ...
                'Press a key to finish!'], ...
                trial_sum,n_trials,perc_corr));
        end

    end

    if lab == false
        Screen('CloseAll');
    end

    ShowCursor(win.h);
    ListenChar(1);

catch e
    ListenChar(1);
    ShowCursor(win.h);
    fclose('all');
    Screen('CloseAll');
    rethrow(e);
    ListenChar(1);

    if use_eyetracker
        % Shutdown Eyelink:
        Eyelink('Shutdown');
    end
end

%% Write results to files using Matlab's amazing new Tables (WELCOME TO THE FUTURE TM):

% sort by trial:
data = sortrows(data, 'trial');
writetable(data, datafilename);

if use_eyetracker
    Eyelink('Command', 'clear_screen 0')
    Eyelink('CloseFile');
    % download data file
    try
        fprintf('Receiving data file ''%s''\n', edf_file_remote );
        status=Eyelink('ReceiveFile',[],eye_data_path,1);
        pause(1);

        movefile([eye_data_path, edf_file_remote], eyedata_fname);
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edf_file_remote, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edf_file_remote, eye_data_path );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', edf_file_remote );
    end
    % Shutdown Eyelink:
    Eyelink('Shutdown');
end

fclose('all');

experiment_duration = toc(experiment_start_time) / 60;
fprintf('The experiment took %2.1f minutes\n', experiment_duration);
fprintf('The observer got %2.1f percent correct in this block', perc_corr);


%% Some subfunctions

    function trial_order=pseudo_rand(data, max_consec)
        % this function adapted from Heiko's `short_trials.m`.
        % max_consec = How many trials between same image allowed >=1 !
        n_t = height(data); % number of trials

        idx = Shuffle(1:n_t);

        accepted = false;

        while ~accepted
            im_vec = data.im_code;
            % FUCK matlab cell arrays...
            im_vec = im_vec(idx);
            for i = 1 : (n_t - max_consec)
                this_im = im_vec(i);

                for j = 1 : max_consec
                    % if the next trial is the same image, jump to a random spot.
                    if  strcmp(im_vec(i + j), this_im)
                        who_jumps = i + j;
                        where_to_jump = randi(n_t);
                        % swap indices:
                        displaced_idx = idx(where_to_jump);
                        idx(where_to_jump) = idx(who_jumps);
                        idx(who_jumps) = displaced_idx;  % swapped with random int.
                    end
                end
            end

            % do a check loop:
            problem = false;
            for i = 1 : (n_t - max_consec)
                this_im = im_vec(i);
                if any(strcmp(im_vec((i+1):(i+max_consec)), this_im))
                    problem = true;
                end
            end

            if ~problem
                accepted = true;
                disp('success in pseudorandom trial ordering');
            else
                disp('trial ordering failed, looping again');
            end
        end
        trial_order = idx;
    end


    function fname = get_natural_im(im_code)
        % return full file path to the image:
        fname =  char(fullfile(orig_im_path, sprintf('%s.png',im_code)));
    end


    function fname = get_model_im(im_code, image_model, version)
        if strcmp(image_model, 'powerspec')
            fname =  char(fullfile(powerspec_path, sprintf(...
                '%s_v%d_syn.png', im_code, version)));
        elseif strcmp(image_model, 'conv5')
            fname =  char(fullfile(conv5_path, sprintf(...
                '%s_version%d.png', im_code, version)));
        elseif strcmp(image_model, 'PS')
            fname =  char(fullfile(ps_path, sprintf('%s_synth_%d.png',im_code,version)));
        end
    end

    function rect = place_rect2(fix_position, radius, angle, rect_base, patch_size)
        % fix position in screen coords, radius (distance from fixation
        % to inner patch edge) in pixels, angle in radians
        % (ccw from right), patch size is diameter of the patch in pixels.
        centre=[(cos(angle) * radius) + fix_position(1), ...
            (sin(angle) * radius) + fix_position(2)+patch_size/2];
        rect = CenterRectOnPoint(rect_base, centre(1), centre(2));
    end

    function rect = place_rect3(fix_position, radius, angle, rect_base, patch_size)
        % fix position in screen coords, radius (distance from fixation
        % to inner patch edge) in pixels, angle in radians
        % (ccw from right), patch size is diameter of the patch in pixels.
        centre=[(cos(angle) * radius) + fix_position(1)+patch_size/2, ...
            (sin(angle) * radius) + fix_position(2)-patch_size/2];
        rect = CenterRectOnPoint(rect_base, centre(1), centre(2));
    end

    function rect = place_rect1(fix_position, radius, angle, rect_base, patch_size)
        % fix position in screen coords, radius (distance from fixation
        % to inner patch edge) in pixels, angle in radians
        % (ccw from right), patch size is diameter of the patch in pixels.
        centre=[(cos(angle) * radius) + fix_position(1)-patch_size/2, ...
            (sin(angle) * radius) + fix_position(2)-patch_size/2];
        rect = CenterRectOnPoint(rect_base, centre(1), centre(2));
    end

    function rect = place_rect(fix_position, radius, angle, rect_base, patch_size)
        % fix position in screen coords, radius (distance from fixation
        % to inner patch edge) in pixels, angle in radians
        % (ccw from right), patch size is diameter of the patch in pixels.
        centre = compute_patch_centre(fix_position, radius, angle, patch_size);
        rect = CenterRectOnPoint(rect_base, centre(1), centre(2));
    end


    function patch_centre = compute_patch_centre(fix_position, radius, angle, patch_size)
        % fix position in screen coords, radius (distance from fixation
        % to inner patch edge) in pixels, angle in radians
        % (ccw from right), patch size is diameter of the patch in pixels.
        patch_centre = [(cos(angle) * (radius + patch_size / 2)) + fix_position(1), ...
            (sin(-angle) * (radius + patch_size / 2)) + fix_position(2)];
    end

    function [tex, im] = make_image_texture(fname, crop_x, crop_y)
        im = im2double(imread(fname));
        end_x = crop_x + params.image_size_px - 1;
        end_y = crop_y + params.image_size_px - 1;

        % crop image according to x, y:
        im = im(crop_y : end_y, crop_x : end_x);

        % add alpha channel (spatial window):
        im(:, :, 2) = spatial_window;

        tex = win.make_texture(im);
    end

% function to check fixation position before trial:
    function check_fixation(check_dist, check_time, recal_time, ...
            ppd, eye_used, fix_position, win, el)
        %%%%%% check initial fixation is within 2deg of the spot %%%%%
        % code here adapted from Will Harrison...
        %
        % check_dist: the distance from fixation that is
        % invalid (in deg).
        % check_time: the duration (ms) to collect valid samples.
        % recal_time: the duration (ms) before triggering recalibration
        % ppd: pixels per degree
        % eye_used: which eye to check
        % fix_position: a tuple containing the fixation position.
        % win: the window object.
        % el: the eyelink object.
        %

        checkFix = 0;

        while checkFix == 0

            % check for presence of a new sample update
            if Eyelink('NewFloatSampleAvailable') > 0

                evt = Eyelink('NewestFloatSample');

                if eye_used ~= -1 % do we know which eye to use yet?
                    % if we do, get current gaze position from sample
                    x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                    y = evt.gy(eye_used+1);

                    fixCheckX = abs(x-(fix_position(1)));
                    fixCheckY = abs(y-fix_position(2));

                    % do we have valid data and is the pupil visible?
                    if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0

                        fixedTimer = GetSecs;
                        tic;
                        while ((fixCheckX < ppd*check_dist) && (fixCheckY < ppd*check_dist));

                            win.draw_fixation(ppd, fix_position, fix_colour_oval, fix_colour_normal);
                            win.flip();

                            evt = Eyelink('NewestFloatSample');

                            x = evt.gx(eye_used+1);
                            y = evt.gy(eye_used+1);

                            fixCheckX = abs(x-(fix_position(1)));
                            fixCheckY = abs(y-fix_position(2));
                            fixTime = GetSecs - fixedTimer;

                            if fixTime > (check_time / 1000)

                                checkFix = 1;
                                break;

                            end

                        end;

                        checkFixTimer = GetSecs;

                        if checkFix < 1

                            while ((fixCheckX > ppd*check_dist) || (fixCheckY > ppd*check_dist));

                                win.draw_fixation(ppd, fix_position, fix_colour_oval, fix_colour_normal);
                                win.flip();

                                evt = Eyelink('NewestFloatSample');
                                x = evt.gx(eye_used+1);
                                y = evt.gy(eye_used+1);
                                fixCheckX = abs(x-(fix_position(1)));
                                fixCheckY = abs(y-fix_position(2));

                                totalCheckFixTime = GetSecs - checkFixTimer;

                                if totalCheckFixTime > (recal_time / 1000)

                                    win.draw_text(['The tracker thinks you are '...
                                        'not looking at the fixation spot. \n\n Recalibrating...']);
                                    win.flip();
                                    WaitSecs(1);

                                    EyelinkDoTrackerSetup(el);

                                    % start recording eye position again
                                    Eyelink('StartRecording');
                                    tic;
                                    checkFixTimer = GetSecs;

                                end;

                            end

                        end

                    end
                else % if we don't, first find eye that's being tracked
                    if 0
                        eye_used = el.RIGHT_EYE;
                    else
                        eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                        if eye_used == el.BINOCULAR; % if both eyes are tracked
                            eye_used = el.LEFT_EYE; % use left eye
                        end
                    end
                end
            end
        end
    end

    function [gaze_pos, eye_used] = return_eye_position(eye_used, el)
        % function to return the current gaze position, checking
        % validity.
        % Note that empty matrix returned for gaze pos if invalid samples.
        % eye_used: which eye to check
        % el: the eyelink object.

        % this series copied from EyelinkExample.m
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            evt = Eyelink( 'NewestFloatSample');
            if eye_used ~= -1 % do we know which eye to use yet?
                % if we do, get current gaze position from sample
                x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used+1);
                % do we have valid data and is the pupil visible?
                if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                    % if data is valid, return gaze position
                    gaze_pos = [x, y];
                else
                    % if data is invalid (e.g. during a blink), return
                    % empty:
                    gaze_pos = [];
                end
            else % if we don't, first find eye that's being tracked
                eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                if eye_used == el.BINOCULAR; % if both eyes are tracked
                    eye_used = el.LEFT_EYE; % use left eye
                end
                gaze_pos = [];
            end
        else
            gaze_pos = [];
        end % if sample available
    end


    function [cum_gaze_pos, eye_used] = cumulative_gaze_pos(cum_gaze_pos, eye_used, el)
        % append gaze pos with newest gaze pos.
        % gaze_pos: an (n_samples, [x, y]) matrix.
        [gaze_pos, eye_used] = return_eye_position(eye_used, el);
        cum_gaze_pos = [cum_gaze_pos; gaze_pos];
    end


    function valid_vec = check_gaze_valid(cum_gaze_pos, fix_position, ok_range_px)
        % check that the gaze positions are within the allowable range.
        % if return value is 1, it's a valid eye sample (0 means gaze
        % position was outside allowable distance).
        dist = sqrt(sum((cum_gaze_pos - repmat(fix_position, size(cum_gaze_pos, 1), 1)).^2, 2));
        valid_vec = ones(size(dist));
        valid_vec(dist > ok_range_px) = 0;
    end


    function [eye_valid, prop_valid] = summarise_trial_eye_valid(valid_vec)
        % returns "eye_valid" == 1 if fixation kept for whole time.
        % returns proportion of valid samples (higher is better).
        eye_valid = all(valid_vec == 1);
        prop_valid = mean(valid_vec);
    end
end
