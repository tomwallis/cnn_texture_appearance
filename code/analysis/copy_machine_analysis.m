
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
crop = crop_size : image_size;

min_offset_orig = 32;
idx_orig = zeros(image_size + crop_size - 1);
crop_orig = crop_size - 32 : image_size + 32;
idx_orig(crop_orig, crop_orig) = 1;
ctr = (image_size + crop_size) / 2 + (-min_offset_orig : min_offset_orig);
idx_orig(ctr, ctr) = 0;

vec = @(x) reshape(x, [], 1);

C = zeros(num_files, num_synth, num_crops+1, num_models);
Corig = zeros(num_files, 1);

for ifile = 1 : num_files
    [~, file_name, ~] = fileparts(files(ifile).name);
    fprintf('Processing %s (%d of %d)\n', file_name, ifile, num_files)

    % load original
    orig = round(double(imread(fullfile(originals_folder, files(ifile).name))) / 256);
    
    Ci = crosscorr(orig, orig(c+1 : end-c, c+1 : end-c));
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
            Ci = crosscorr(orig, synth(c+1 : end-c, c+1 : end-c));
            C(ifile, isynth, 1, imodel) = max(vec(Ci(crop, crop)));
            for icrop = 1 : num_crops
                ij = fix(rand(1, 2) * (image_size - crop_size));
                Ci = crosscorr(orig, synth(ij(1)+1 : ij(1)+crop_size, ij(2)+1 : ij(2)+crop_size));
                C(ifile, isynth, 1+icrop, imodel) = max(vec(Ci(crop, crop)));
            end
        end
    end
end
fprintf('\n')

%% write csv
fid = fopen(fullfile(originals_folder, 'data.csv'), 'w');
fprintf(fid, repmat('%s, ', 1, num_models), models{:});
fprintf(fid, 'original\n');
for i = 1 : num_files
    for j = 1 : num_models
        fprintf(fid, '%.8f, ', Cm(i, j));
    end
    fprintf(fid, '%.8f\n', Corig(i));
end
fclose(fid);
