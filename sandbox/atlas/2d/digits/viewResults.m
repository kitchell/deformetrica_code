close all;

scale = 1;
indexes = 1:5;

% load the initial & final templates
Itemp_i = imread('data/digit_2_mean.png');
Itemp_f = imread('output/BayesianAtlas_digit_2_mean.png');

% % % Itemp_f_sr = importdata('../Output/Atlas_digit_2_mean.txt');
% % % Itemp_f_sr = round((Itemp_f_sr - min(Itemp_f_sr(:))) * 255);

a = min(Itemp_i(:)); b = max(Itemp_i(:));
Itemp_i = (Itemp_i-a)*(255/double(b-a));
Itemp_f = (Itemp_f-a)*(255/double(b-a));

% % % Itemp_i = (Itemp_i-min(Itemp_i(:)))*(255/double(max(Itemp_i(:))-min(Itemp_i(:))));
% % % Itemp_f = (Itemp_f-min(Itemp_f(:)))*(255/double(max(Itemp_f(:))-min(Itemp_f(:))));

% load target data and results
Itarget = cell(1,5);
Iresult = cell(1,5);
for s = 1:5
	Itarget{s} = imread(['data/digit_2_sample_', num2str(indexes(s)), '.png']);
    Iresult{s} = imread(['output/BayesianAtlas_digit_2_mean_to_subject_', num2str(indexes(s)-1), '__t_9.png']);

    Iresult{s} = (Iresult{s}-a)*(255/double(b-a));
    
% % %     Iresult{s} = (Iresult{s}-min(Iresult{s}(:)))*(255/double(max(Iresult{s}(:))-min(Iresult{s}(:))));
% % %     Itarget{s} = (Itarget{s}-min(Itarget{s}(:)))*(255/double(max(Itarget{s}(:))-min(Itarget{s}(:))));
end

% plot
figure;
subplot(2,6,1)
imshow(imresize(Itemp_i, scale));
title('Initial template');
subplot(2,6,7)
imshow(imresize(Itemp_f, scale));
%imshow(Itemp_f_sr, [0, 255]);
title('Final template');
for s=1:5
	subplot(2,6,1+s);
    imshow(imresize(Itarget{s}, scale));
    title(['Target ', num2str(indexes(s))]);
    
    subplot(2,6,7+s);
    imshow(imresize(Iresult{s}, scale));
    title(['Result ', num2str(indexes(s))]);
end
set(gcf,'OuterPosition',[-1500 1750 1750 1000]);


% % % % plots the physical coordinates
% % % [xpoff, ypoff] = textread('../Output/0_OFF_Y0.txt', '%f %f');
% % % [xpon, ypon] = textread('../Output/0_ON_Y0.txt', '%f %f');
% % % 
% % % figure; 
% % % subplot(2,1,1);
% % % plot(1:length(xoff), xoff, 'r'); hold on; 
% % % plot(1:length(xon), xon, 'b');
% % % xlim([1, 255]);
% % % title('Physical x coordinates. Parametric = blue ; LinearInterp = red.');
% % % 
% % % subplot(2,1,2);
% % % plot(0:255, yoff, 'r'); hold on; 
% % % plot(0:255, yon, 'b');
% % % xlim([0, 255]);
% % % title('Physical y coordinates. Parametric = blue ; LinearInterp = red.');
% % % set(gcf,'OuterPosition',[-1500 1750 1750 1000]);
% % % 
% % % % plots the spatial gradient of the image
% % % [xoff, yoff] = textread('../Output/0_OFF_NablaI.txt', '%f %f');
% % % [xon, yon] = textread('../Output/0_ON_NablaI.txt', '%f %f');
% % % 
% % % figure; 
% % % subplot(2,1,1);
% % % plot((0:(length(xoff)-1))/length(xoff), xoff, 'r'); hold on; 
% % % plot((0:(length(xon)-1))/length(xon), xon, 'b');
% % % title('Image gradient on x. Parametric = blue ; LinearInterp = red.');
% % % 
% % % subplot(2,1,2);
% % % plot((0:(length(yoff)-1))/length(yoff), yoff, 'r'); hold on; 
% % % plot((0:(length(yon)-1))/length(yon), yon, 'b');
% % % title('Image gradient on y. Parametric = blue ; LinearInterp = red.');
% % % set(gcf,'OuterPosition',[-1500 1750 1750 1000]);

% % % % plots the spatial gradient of the match
% % % [xmoff, ymoff] = textread('../Output/0_OFF_GradMatch.txt', '%f %f');
% % % [xmon, ymon] = textread('../Output/0_ON_GradMatch.txt', '%f %f');
% % % 
% % % figure; 
% % % subplot(2,1,1);
% % % plot((0:(length(xmoff)-1))/length(xmoff), xmoff, 'r'); hold on; 
% % % plot((0:(length(xmon)-1))/length(xmon), xmon, 'b');
% % % title('Match gradient on x. Parametric = blue ; LinearInterp = red.');
% % % 
% % % subplot(2,1,2);
% % % plot((0:(length(ymoff)-1))/length(ymoff), ymoff, 'r'); hold on; 
% % % plot((0:(length(ymon)-1))/length(ymon), ymon, 'b');
% % % title('Match gradient on y. Parametric = blue ; LinearInterp = red.');
% % % set(gcf,'OuterPosition',[-1500 1750 1750 1000]);

% % % % plots the intensity gradient
% % % [ioff] = textread('../Output/0_OFF_NablaIntensity.txt', '%f');
% % % [ion] = textread('../Output/0_ON_NablaIntensity.txt', '%f');
% % % 
% % % figure; 
% % % plot((0:(length(ioff)-1))/length(ioff), ioff, 'r'); hold on; 
% % % plot((0:(length(ion)-1))/length(ion), ion, 'b');
% % % title('Intensity gradient. Parametric = blue ; LinearInterp = red.');
% % % set(gcf,'OuterPosition',[-1500 1750 1750 1000]);
% % % 
% % % % plots the CP gradient 
% % % [xcpoff, ycpoff] = textread('../Output/0_OFF_NablaCP.txt', '%f %f');
% % % [xcpon, ycpon] = textread('../Output/0_ON_NablaCP.txt', '%f %f');
% % % 
% % % figure; 
% % % subplot(2,1,1);
% % % plot((0:(length(xcpoff)-1))/length(xcpoff), xcpoff, 'r'); hold on; 
% % % plot((0:(length(xcpon)-1))/length(xcpon), xcpon, 'b');
% % % title('CP gradient on x. Parametric = blue ; LinearInterp = red.');
% % % 
% % % subplot(2,1,2);
% % % plot((0:(length(ycpoff)-1))/length(ycpoff), ycpoff, 'r'); hold on; 
% % % plot((0:(length(ycpon)-1))/length(ycpon), ycpon, 'b');
% % % title('CP gradient on y. Parametric = blue ; LinearInterp = red.');
% % % set(gcf,'OuterPosition',[-1500 1750 1750 1000]);
% % % 
% % % % plots the momenta gradient
% % % [xmomoff, ymomoff] = textread('../Output/0_OFF_NablaMom.txt', '%f %f');
% % % [xmomon, ymomon] = textread('../Output/0_ON_NablaMom.txt', '%f %f');
% % % 
% % % figure; 
% % % subplot(2,1,1);
% % % plot((0:(length(xmomoff)-1))/length(xmomoff), xmomoff, 'r'); hold on; 
% % % plot((0:(length(xmomon)-1))/length(xmomon), xmomon, 'b');
% % % title('Mom gradient on x. Parametric = blue ; LinearInterp = red.');
% % % 
% % % subplot(2,1,2);
% % % plot((0:(length(ymomoff)-1))/length(ymomoff), ymomoff, 'r'); hold on; 
% % % plot((0:(length(ymomon)-1))/length(ymomon), ymomon, 'b');
% % % title('Mom gradient on y. Parametric = blue ; LinearInterp = red.');
% % % set(gcf,'OuterPosition',[-1500 1750 1750 1000]);





