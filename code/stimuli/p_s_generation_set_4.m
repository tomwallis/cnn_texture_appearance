% Matlab script to run texture synthesization on patches.
%
% tsawallis wrote it.

% add the pyramid and texture toolboxes to matlab's search path:
% addpath(genpath('/home/tomw/matlab_toolboxes/matlabPyrTools'))
% addpath(genpath('/home/tomw/matlab_toolboxes/matlabPyrTools/MEX'))
% addpath(genpath('/home/tomw/matlab_toolboxes/textureSynth'))
% addpath(genpath('/home/tomw/matlab_toolboxes/textureSynth/MEX'))

s = RandStream('mt19937ar', 'Seed', 557734);
RandStream.setGlobalStream(s);

% paths to images:

% (because matlab sucks at strings and directories, this code sucks:
this_dir = pwd;
top_dir = this_dir(1:end-12);
raw_path = [top_dir, 'stimuli/stimulus_set_4/preprocessed_ims'];
out_path = [top_dir, 'stimuli/stimulus_set_4/ps_synths'];

if ~exist(out_path, 'dir')
    mkdir(out_path);
end

% ims:
ims = dir([raw_path, '/*.png']);

%% Parameters

n_unique = 10; % the number of unique textures to save for each original patch.

Nor = 4; % Number of orientations
Nsc = 4; % max number of scales for 256 images = 6.
Na = 11;  % Spatial neighborhood is Na x Na coefficients
% It must be an odd number!

Niter = 50;	% Number of iterations of synthesis loop

n_par_cores = 4;  % the number of compute cores available.

%% normal for loop

% for i = 1 : length(ims)
% for i = 1 : 1
%     file = ims(i).name;
%     orig_im = imread([raw_path, '/', file]);
%     im = im2double(orig_im);
%     
%     im_size = size(im, 1);
%     
%     Nsx = size(im, 1);	% Size of synthetic image is Nsy x Nsx
%     Nsy = size(im, 1);	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)
%     
%     params = textureAnalysis(im, Nsc, Nor, Na);
%     
%     for j = 1 : n_unique
%         success = 0;
%         while success == 0
%             try
%                 res = textureSynthesis(params, [Nsy Nsx], Niter);
%                 
%                 fname = sprintf('%s/%s_synth_%s.png', out_path, ...
%                     ims(i).name(1:end-4), num2str(j));
%                 
%                 imwrite(res, fname, 'bitdepth', 16);
%                 
%                 success = 1;
%                 
%             catch res_failed
%                 
%             end
%         end
%     end
% end

%% Parallel for loop, if required:

try
    % open parallel pool:
    pool = parpool(n_par_cores);
    
    parfor i = 1 : length(ims)
        file = ims(i).name;
        orig_im = imread([raw_path, '/', file]);
        im = im2double(orig_im);
               
        im_size = size(im, 1);
        
        Nsx = size(im, 1);	% Size of synthetic image is Nsy x Nsx
        Nsy = size(im, 1);	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)
        
        params = textureAnalysis(im, Nsc, Nor, Na);
        
        for j = 1 : n_unique
            success = 0;
            while success == 0
                try
                    res = textureSynthesis(params, [Nsy Nsx], Niter);
                    
                    fname = sprintf('%s/%s_synth_%s.png', out_path, ...
                        ims(i).name(1:end-4), num2str(j));
                    
                    % save as a 16 bit image:
                    imwrite(res, fname, 'bitdepth', 16);
                    
                    success = 1;
                    
                catch res_failed
                    
                end
            end
        end
    end
    
    delete(gcp('nocreate'))
    %     matlabpool close
catch ME
    delete(gcp('nocreate'))
    %     matlabpool close
    
end


