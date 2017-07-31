% Alex Ecker wrote it. Tom modified some file paths.

[folder, ~] = fileparts(mfilename('fullpath'));
image_folder = fullfile(folder, '../../stimuli/stimulus_set_4');
originals_folder = fullfile(image_folder, 'preprocessed_ims');
output_folder = fullfile(folder, '../../results/copy_machine_crosscor');

files = dir(fullfile(originals_folder, '*.png'));
num_files = numel(files);

models = {'conv5', 'powerspec', 'ps_synths'};
num_models = numel(models);

rng(1);

num_synth = 10;
image_size = 256;
crop_size = 128;
c = crop_size / 2;
num_crops = 10;

min_offset_orig = 32;
idx_orig = ones(image_size - crop_size + 1);
ctr = (image_size - crop_size) / 2 + 1 + (-min_offset_orig : min_offset_orig);
idx_orig(ctr, ctr) = 0;

C = zeros(num_files, num_synth, num_crops+1, num_models);
Corig = zeros(num_files, 1);

for ifile = 1 : num_files
    [~, file_name, ~] = fileparts(files(ifile).name);
    fprintf('Processing %s (%d of %d)\n', file_name, ifile, num_files)

    % load original
    orig = round(double(imread(fullfile(originals_folder, files(ifile).name))) / 256);
    
    ctr = image_size / 2;
    orig_center_crop = orig(ctr-c+1 : ctr+c, ctr-c+1 : ctr+c);
    Ci = crosscorrcoef(orig, orig_center_crop);
    Corig(ifile) = max(Ci(idx_orig > 0));
    
    for imodel = 1 : num_models
        fprintf('  %s\n', models{imodel})

        % load synthesized images
        synth_folder = fullfile(image_folder, models{imodel});
        synth_files = dir(fullfile(synth_folder, [file_name '*']));
        for isynth = 1 : num_synth
            synth = double(imread(fullfile(synth_folder, synth_files(isynth).name)));
            if mean(synth(:)) > 256
                synth = round(synth / 256);
            end
            synth_center_crop = synth(ctr-c+1 : ctr+c, ctr-c+1 : ctr+c);
            Ci = crosscorrcoef(orig, synth_center_crop);
            C(ifile, isynth, 1, imodel) = max(Ci(:));
            for icrop = 1 : num_crops
                ij = fix(rand(1, 2) * (image_size - crop_size));
                synth_rand_crop = synth(ij(1)+1 : ij(1)+crop_size, ij(2)+1 : ij(2)+crop_size);
                Ci = crosscorrcoef(orig, synth_rand_crop);
                C(ifile, isynth, 1+icrop, imodel) = max(Ci(:));
            end
        end
    end
end
fprintf('\n')

%% write csv
Cm = squeeze(mean(mean(C, 2), 3));
fid = fopen(fullfile(originals_folder, sprintf('data%d.csv', crop_size)), 'w');
fprintf(fid, repmat('%s, ', 1, num_models), models{:});
fprintf(fid, 'original\n');
for i = 1 : num_files
    for j = 1 : num_models
        fprintf(fid, '%.8f, ', Cm(i, j));
    end
    fprintf(fid, '%.8f\n', Corig(i));
end
fclose(fid);
