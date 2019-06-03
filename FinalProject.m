clear all; close all;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Part A: Synthesis
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%A.1 Generating Gaussian Noise from Uniform Noise
Nx = 128;
My = 128;
N = Nx * My;
M = 50;

%randomSeeds = randi(500,1,6);
randomseeds = [27,89,193,460,495,239];

img = {};
U = {};

%A.1.1
for i = 1:6
    currrng = rng(randomseeds(i))
    img{i} = rand(Nx, My);
    figure(i);
    imagesc(img{i});
    colormap gray
    axis image;
    U{i} = reshape(img{i}, [N,1]);
end

figure(30);
h = [];
title("U Distributions");
hold on;
h(1) = histogram(reshape(U{1}, [N, 1]), 10, 'DisplayName','u1');
h(2) = histogram(reshape(U{2}, [N, 1]), 10, 'DisplayName','u2');
h(3) = histogram(reshape(U{3}, [N, 1]), 10, 'DisplayName','u3');
h(4) = histogram(reshape(U{4}, [N, 1]), 10, 'DisplayName','u4');
h(5) = histogram(reshape(U{5}, [N, 1]), 10, 'DisplayName','u5');
hold off;
legend(h);

%A.1.2
%SEE ATTACHED DOCUMENT SECTION 2.1

%A.1.3
g = {};
for i = 1:6
    g{i} = zeros([N, 1]);
    for n = 0:(N/2-1)
        g{i}(2*n+1) = sqrt(-2*log(U{i}(2*n+1)))*cos(2*pi*U{i}(2*n+2));
        g{i}(2*n+2) = sqrt(-2*log(U{i}(2*n+1)))*sin(2*pi*U{i}(2*n+2));
    end
    stats = [mean(g{i}) std(g{i}) var(g{i})]
end

figure(31);
h = [];
title("g Distributions");
hold on;
h(1) = histogram(reshape(g{1}, [N, 1]), 10, 'DisplayName','g1');
h(2) = histogram(reshape(g{2}, [N, 1]), 10, 'DisplayName','g2');
h(3) = histogram(reshape(g{3}, [N, 1]), 10, 'DisplayName','g3');
h(4) = histogram(reshape(g{4}, [N, 1]), 10, 'DisplayName','g4');
h(5) = histogram(reshape(g{5}, [N, 1]), 10, 'DisplayName','g5');
hold off;
legend(h);

%A.1.4
s = {};
for i = 1:5
    s{i} = zeros([N, 1]);
    for n = 1:i+1
        s{i} = s{i} + U{n}
    end

    s{i} = reshape(s{i}, [Nx,My]);
%     figure(i);
%     imagesc(s{i});
%     colormap gray
%     axis image;
end

%A.2 Generating Control Noise from Deterministic Data
%A.2.1
targetSelectionBmp = double(imread('Target/target0.bmp'));
Ar=targetSelectionBmp(:,:,1);
Ag=targetSelectionBmp(:,:,2);
Ab=targetSelectionBmp(:,:,3);
targetSelection=(Ar+Ag+Ab)/3;
clutterSelectionBmp = double(imread('Clutter/Clutter0.bmp'));
Ar=clutterSelectionBmp(:,:,1);
Ag=clutterSelectionBmp(:,:,2);
Ab=clutterSelectionBmp(:,:,3);
clutterSelection=(Ar+Ag+Ab)/3;

Ht = abs(fft2(targetSelection));
Hc = abs(fft2(clutterSelection));

figure(1);
imagesc(abs(ifftshift(Ht)));
title('Shifted fft of Target Image');
colormap gray
axis image;

figure(2);
imagesc(abs(ifftshift(Hc)));
title('Shifted fft of Clutter Image');
colormap gray
axis image;


%A.2.2
At = 0 + (2*pi-0).*rand(N,1);
Ac = 0 + (2*pi-0).*rand(N,1);

Ht = reshape(Ht, [N,1]);
Hc = reshape(Hc, [N,1]);

for i = 1:N
    At(i) = Ht(i)*exp(j*At(i));
    Ac(i) = Hc(i)*exp(j*Ac(i));
end

At = reshape(At, [Nx,My]);
Ac = reshape(Ac, [Nx,My]);
Ht = reshape(Ht, [Nx,My]);
Hc = reshape(Hc, [Nx,My]);

%apply symmetry rules
for n = 1:My
    for m = 1:Nx
        if (m == 1 && n > My/2)
            n1 = 2+My-n;            
            At(m,n) = real(At(m,n1))-1i*imag(At(m,n1));
            Ac(m,n) = real(Ac(m,n1))-1i*imag(Ac(m,n1));
        end
        if (n == 1 && m > Nx/2)
            m1 = 2+Nx-m;
            At(m,n) = real(At(m1,n))-1i*imag(At(m1,n));
            Ac(m,n) = real(Ac(m1,n))-1i*imag(Ac(m1,n));
        end
    end
end

imageTarget2 = ifft2(At, 'symmetric');
imageClutter2 = ifft2(Ac, 'symmetric');

figure(3);
imagesc(imageTarget2);
title('Uniform Control Noise From Target Image');
colormap gray
axis image;

figure(4);
imagesc(imageClutter2);
title('Uniform Control Noise From Clutter Image');
colormap gray
axis image;


%A.2.3
gt = normrnd(0, 1, [N,1]);
stats = [mean(gt) std(gt) var(gt)]
gc = normrnd(0, 1, [N,1]);
stats = [mean(gc) std(gc) var(gc)]
gt = reshape(gt, [Nx,My]);
gc = reshape(gc, [Nx,My]);

Gt = fft2(gt).*Ht;
Gc = fft2(gc).*Hc;

%apply symmetry rules
for n = 1:My
    for m = 1:Nx
        if (m == 1 && n > My/2)
            n1 = 2+My-n;            
            Gt(m,n) = real(Gt(m,n1))-1i*imag(Gt(m,n1));
            Gc(m,n) = real(Gc(m,n1))-1i*imag(Gc(m,n1));
        end
        if (n == 1 && m > Nx/2)
            m1 = 2+Nx-m;
            Gt(m,n) = real(Gt(m1,n))-1i*imag(Gt(m1,n));
            Gc(m,n) = real(Gc(m1,n))-1i*imag(Gc(m1,n));
        end
    end
end

imageTarget3 = ifft2(Gt, 'symmetric');
imageClutter3 = ifft2(Gc, 'symmetric');

figure(5);
imagesc(imageTarget3);
title('Normal Control Noise From Target Image');
colormap gray
axis image;

figure(6);
imagesc(imageClutter3);
title('Normal Control Noise From Clutter Image');
colormap gray
axis image;

%A.2.4
imageTarget2PSD = log10(abs(fftshift(fft2(double(255*mat2gray(imageTarget2))))).^2);
imageClutter2PSD = log10(abs(fftshift(fft2(double(255*mat2gray(imageClutter2))))).^2);
imageTarget3PSD = log10(abs(fftshift(fft2(double(255*mat2gray(imageTarget3))))).^2);
imageClutter3PSD = log10(abs(fftshift(fft2(double(255*mat2gray(imageClutter3))))).^2);
HtPSD = log10(abs(fftshift(fft2(double(255*mat2gray(targetSelection))))).^2);
HcPSD = log10(abs(fftshift(fft2(double(255*mat2gray(clutterSelection))))).^2);

imageTarget2AutoCorr = normxcorr2(double(255*mat2gray(imageTarget2)), double(255*mat2gray(imageTarget2)));
imageClutter2AutoCorr = normxcorr2(double(255*mat2gray(imageClutter2)), double(255*mat2gray(imageClutter2)));
imageTarget3AutoCorr = normxcorr2(double(255*mat2gray(imageTarget3)), double(255*mat2gray(imageTarget3)));
imageClutter3AutoCorr = normxcorr2(double(255*mat2gray(imageClutter3)), double(255*mat2gray(imageClutter3)));
HtAutoCorr = normxcorr2(double(255*mat2gray(targetSelection)), double(255*mat2gray(targetSelection)));
HcAutoCorr = normxcorr2(double(255*mat2gray(clutterSelection)), double(255*mat2gray(clutterSelection)));


figure(7);
imagesc(imageTarget2PSD);
title('Filtered At PSD');
colormap gray
axis image;

figure(8);
imagesc(imageClutter2PSD);
title('Filtered Ac PSD');
colormap gray
axis image;

figure(9);
imagesc(imageTarget3PSD);
title('Filtered gt PSD');
colormap gray
axis image;

figure(10);
imagesc(imageClutter3PSD);
title('Filtered gc PSD');
colormap gray
axis image;

figure(17);
imagesc(HtPSD);
title('Ht PSD');
colormap gray
axis image;

figure(18);
imagesc(HcPSD);
title('Hc PSD');
colormap gray
axis image;

immse(imageTarget3PSD, imageTarget2PSD)
immse(imageClutter3PSD, imageClutter2PSD)
immse(imageTarget3PSD, HtPSD)
immse(imageClutter3PSD, HcPSD)
immse(HtPSD, imageTarget2PSD)
immse(HcPSD, imageClutter2PSD)

immse(imageTarget3AutoCorr, imageTarget2AutoCorr)
immse(imageClutter3AutoCorr, imageClutter2AutoCorr)
immse(imageTarget3AutoCorr, HtAutoCorr)
immse(imageClutter3AutoCorr, HcAutoCorr)
immse(HtAutoCorr, imageTarget2AutoCorr)
immse(HcAutoCorr, imageClutter2AutoCorr)
%SEE ATTACHED DOCUMENT 




%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Part B: Analysis
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%B.1.a
tI = {};
tI{1} = g{1} + 2;
tI{2} = g{2} + 3;
tI{3} = g{3} + 4;

cI = {};
cI{1} = g{4} + 1;
cI{2} = g{5} + 1;
cI{3} = g{6} + 1;

%B.1.b
t ={};
t{1} = 4.*g{1} + g{2} + 2.*g{3} + 28;
t{2} = 2.*g{1} + 4.*g{2} + g{3} + 10;
t{3} = 2.*g{1} + g{2} + 4.*g{3} + 21;

c = {};
c{1} = 4.*g{4} + 2.*g{5} + g{6} + 7;
c{2} = g{4} + 4.*g{5} + 2.*g{6} + 7;
c{3} = 2.*g{4} + g{5} + 4.*g{6} + 7;

%B.2

figure(11);
h = [];
title("PlotB - 1 - Gray Image Intensity");
hold on;
H = hist(double(255*mat2gray(U{1})), M);
h(1) = plot(H, 'DisplayName','u1');
H = hist(double(255*mat2gray(g{1})), M);
h(2) = plot(H, 'DisplayName','g1');
hold off;
legend(h);

figure(32);
h = [];
title("PlotB - 1 - Value Histogram");
hold on;
h(1) = histogram(U{1}, M, 'DisplayName','u1');
h(2) = histogram(g{1}, M, 'DisplayName','g1');
hold off;
legend(h);

figure(12);
h = [];
title("PlotB - 2 - Gray Image Intensity");
hold on;
H = hist(double(255*mat2gray(reshape(s{1}, [N, 1]))), M);
h(1) = plot(H, 'DisplayName','s1');
H = hist(double(255*mat2gray(reshape(s{2}, [N, 1]))), M);
h(2) = plot(H, 'DisplayName','s2');
H = hist(double(255*mat2gray(reshape(s{3}, [N, 1]))), M);
h(3) = plot(H, 'DisplayName','s3');
H = hist(double(255*mat2gray(reshape(s{4}, [N, 1]))), M);
h(4) = plot(H, 'DisplayName','s4');
H = hist(double(255*mat2gray(reshape(s{5}, [N, 1]))), M);
h(5) = plot(H, 'DisplayName','s5');
hold off;
legend(h);

figure(33);
h = [];
title("PlotB - 2 - Value Histogram");
hold on;
h(1) = histogram(s{1}, M, 'DisplayName','s1');
h(2) = histogram(s{2}, M, 'DisplayName','s2');
h(3) = histogram(s{3}, M, 'DisplayName','s3');
h(4) = histogram(s{4}, M, 'DisplayName','s4');
h(5) = histogram(s{5}, M, 'DisplayName','s5');
hold off;
legend(h);

figure(13);
h = [];
title("PlotB - 3 - Gray Image Intensity");
hold on;
H = hist(double(255*mat2gray(tI{1})), M);
h(1) = plot(H, 'DisplayName','tI1');
H = hist(double(255*mat2gray(cI{1})), M);
h(2) = plot(H, 'DisplayName','cI1');
hold off;
legend(h);

figure(25);
h = [];
title("PlotB - 3 - Value Histogram");
hold on;
h(1) =  histogram(tI{1}, M, 'DisplayName','tI1');
h(2) =  histogram(cI{1}, M, 'DisplayName','cI1');
hold off;
legend(h);

figure(14);
h = [];
title("PlotB - 4 - Gray Iamge Intensity");
hold on;
H = hist(double(255*mat2gray(tI{2})), M);
h(1) = plot(H, 'DisplayName','tI2');
H = hist(double(255*mat2gray(cI{2})), M);
h(2) = plot(H, 'DisplayName','cI2');
hold off;
legend(h);

figure(26);
h = [];
title("PlotB - 4 - Value Histogram");
hold on;
h(1) =  histogram(tI{2}, M, 'DisplayName','tI2');
h(2) =  histogram(cI{2}, M, 'DisplayName','cI2');
hold off;
legend(h);

figure(15);
h = [];
title("PlotB - 5 - Gray Iamge Intensity");
hold on;
H = hist(double(255*mat2gray(tI{3})), M);
h(1) = plot(H, 'DisplayName','tI3');
H = hist(double(255*mat2gray(cI{3})), M);
h(2) = plot(H, 'DisplayName','cI3');
hold off;
legend(h);

figure(27);
h = [];
title("PlotB - 5 - Value Histogram");
hold on;
h(1) =  histogram(tI{3}, M, 'DisplayName','tI3');
h(2) =  histogram(cI{3}, M, 'DisplayName','cI3');
hold off;
legend(h);

%B.3
Kt = zeros(3,3);
for i = 1:3
    for j = 1:3
        mu_n(1:N,1) = mean(t{j});
        mu_m(1:N,1) = mean(t{i});
        Kt(i,j) = (1/(N-1))*(t{i}-mu_m)'*(t{j}-mu_n);
    end
end

Kc = zeros(3,3);
for i = 1:3
    for j = 1:3
        mu_n(1:N,1) = mean(c{j});
        mu_m(1:N,1) = mean(c{i});
        Kc(i,j) = (1/(N-1))*(c{i}-mu_m)'*(c{j}-mu_n);
    end
end

KtI = zeros(3,3);
for i = 1:3
    for j = 1:3
        mu_n(1:N,1) = mean(tI{j});
        mu_m(1:N,1) = mean(tI{i});
        KtI(i,j) = (1/(N-1))*(tI{i}-mu_m)'*(tI{j}-mu_n);
    end
end

KcI = zeros(3,3);
for i = 1:3
    for j = 1:3
        mu_n(1:N,1) = mean(cI{j});
        mu_m(1:N,1) = mean(cI{i});
        KcI(i,j) = (1/(N-1))*(cI{i}-mu_m)'*(cI{j}-mu_n);
    end
end


%B.4
mu_t = zeros(3,1);
mu_t(1) = mean(tI{1});
mu_t(2) = mean(tI{2});
mu_t(3) = mean(tI{3});

mu_c = zeros(3,1);
mu_c(1) = mean(cI{1});
mu_c(2) = mean(cI{2});
mu_c(3) = mean(cI{3});


%B.5

tempTargCopy = targetSelection(:);
%flush out zeros so we can get a good look at the rest of the
%image that is not zeros.
nonzeroIndices = tempTargCopy > 0;
tempTargCopy = tempTargCopy(nonzeroIndices);

tempNoiseCopy = imageTarget3(:);
%normalize to grayscale
tempNoiseCopy = double(255*mat2gray(tempNoiseCopy));
%flush out zeros so we can get a good look at the rest of the
%image that is not zeros.
nonzeroIndices = tempNoiseCopy > 0;
tempNoiseCopy = tempNoiseCopy(nonzeroIndices);

%plot the two on the same figure
figure(16);
h = [];
hold on;
%tempH1 = hist(tempTargCopy, M);
tempH1 = histogram(tempTargCopy, M, 'Normalization','probability', 'DisplayName','target0 Intensity');
%tempH2 = hist(tempNoiseCopy, M);
tempH2 = histogram(tempNoiseCopy, M, 'Normalization','probability', 'DisplayName','control noise Intensity');
h(1) = tempH1;
h(2) = tempH2;
%h(1) = plot(tempH1, 'DisplayName','target0 Intensity');
%h(2) = plot(tempH2, 'DisplayName','control noise Intensity');
title("Target and Control Noise Intensity Histograms");
hold off;
legend(h);

%SEE ATTACHED DOCUMENT 




%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Part C: Detection and Discrimination
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%C.1
%find Ji for tIi and cIi
JI = {};
for i = 1:3
    sig_t_squared = 0;
    sig_c_squared = 0;
    
    for j = 1:N
        sig_t_squared = sig_t_squared + (tI{i}(j) - mu_t(i))^2;
        sig_c_squared = sig_c_squared + (cI{i}(j) - mu_c(i))^2;
    end
    sig_t_squared = sig_t_squared*(1/(N-1));
    sig_c_squared = sig_c_squared*(1/(N-1));
    
    JI{i} = (mu_t(i) - mu_c(i))^2/(sig_t_squared + sig_c_squared);
end

%find Ji for ti and ci
mu_t_notI = zeros(3,1);
mu_t_notI(1) = mean(t{1});
mu_t_notI(2) = mean(t{2});
mu_t_notI(3) = mean(t{3});

mu_c_notI = zeros(3,1);
mu_c_notI(1) = mean(c{1});
mu_c_notI(2) = mean(c{2});
mu_c_notI(3) = mean(c{3});
J = {};
for i = 1:3
    sig_t_squared = 0;
    sig_c_squared = 0;
    
    for j = 1:N
        sig_t_squared = sig_t_squared + (t{i}(j) - mu_t_notI(i))^2;
        sig_c_squared = sig_c_squared + (c{i}(j) - mu_c_notI(i))^2;
    end
    sig_t_squared = sig_t_squared*(1/(N-1));
    sig_c_squared = sig_c_squared*(1/(N-1));
    
    J{i} = (mu_t_notI(i) - mu_c_notI(i))^2/(sig_t_squared + sig_c_squared);
end

%C.2
%C.2.a
%forming MLR test with tI{3} and cI{3} as they correspond to the largest J
%value
K_avg = (KtI(3,3)+KcI(3,3))./2;
mu_t_sample = mu_t(3);
mu_c_sample = mu_c(3);
eta = 0.5*(mu_t_sample'*inv(K_avg)*mu_t_sample - mu_c_sample'*inv(K_avg)*mu_c_sample);
I_one_rvs_target_response = [];
I_one_rvs_clutter_response = [];
I_one_rvs_prob_false_alarm = 0;
I_one_rvs_prob_miss = 0;
for i = 1:N
    currTargetCalc = tI{3}(i)'*(inv(K_avg)*(mu_t_sample - mu_c_sample));
    currClutterCalc = cI{3}(i)'*(inv(K_avg)*(mu_t_sample - mu_c_sample));
    I_one_rvs_target_response = horzcat(I_one_rvs_target_response, currTargetCalc);
    I_one_rvs_clutter_response = horzcat(I_one_rvs_clutter_response, currClutterCalc);
    
    if currTargetCalc < eta
        I_one_rvs_prob_miss = I_one_rvs_prob_miss + 1;
    end
    if currClutterCalc > eta
        I_one_rvs_prob_false_alarm = I_one_rvs_prob_false_alarm + 1;
    end
end
I_one_rvs_prob_false_alarm = I_one_rvs_prob_false_alarm / N;
I_one_rvs_prob_miss = I_one_rvs_prob_miss / N;
I_one_rvs_prob_error = I_one_rvs_prob_false_alarm + I_one_rvs_prob_miss;
I_one_rvs_decision_boundary(1,1:N) = eta;

%"dumb" minimization loop to look for MPE(an improvised minimum walking algorithm)
test_eta = [eta];
Prob_errors = [I_one_rvs_prob_error];
for i = 1:5
    degree = 10^-i;
    
    %carry over old candidates
    new_test_eta = test_eta;
    new_Prob_errors = Prob_errors;
    
    for k = 1:length(test_eta)
        currEta = test_eta(k);
        % test back walk
        for j = 1:10
            check_eta = currEta-j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);

            test_back = (sum(I_one_rvs_target_response < check_eta)+sum(I_one_rvs_clutter_response > check_eta))/N;
            new_Prob_errors = horzcat(new_Prob_errors, test_back);
        end

        %test forward walk
        for j = 1:10
            check_eta = currEta+j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);
            test_forward = (sum(I_one_rvs_target_response < check_eta)+sum(I_one_rvs_clutter_response > check_eta))/N; 
            new_Prob_errors = horzcat(new_Prob_errors, test_forward);
        end

        %take our 10 best candidates
        [Prob_errors, errorSort] = sort(new_Prob_errors);
        test_eta = new_test_eta(errorSort);
        if (length(test_eta) > 10)
           Prob_errors = Prob_errors(1:10);
           test_eta = test_eta(1:10);
        end
    end
end
I_one_rvs_min_prob_error = Prob_errors(1);
I_one_rvs_decision_boundary_MPE(1,1:N) = test_eta(1);
sprintf("Probabiliy hunting for one independent r.v.s. found a MPE of %f at a boundary of %f compared to the orginal guess eta of %f's PE of %f",I_one_rvs_min_prob_error, test_eta(1), eta, I_one_rvs_prob_error)



%forming MLR test with tI{2}, tI{3}, cI{2}, and cI{3} as they correspond to the two largest J
%values
K_avg = (KtI(2:3,2:3)+KcI(2:3,2:3))./2;
mu_t_sample = mu_t(2:3);
mu_c_sample = mu_c(2:3);
eta = 0.5*(mu_t_sample'*inv(K_avg)*mu_t_sample - mu_c_sample'*inv(K_avg)*mu_c_sample);
I_two_rvs_target_response = [];
I_two_rvs_clutter_response = [];
t_sample = [tI{2} tI{3}]';
c_sample = [cI{2} cI{3}]';
I_two_rvs_prob_false_alarm = 0;
I_two_rvs_prob_miss = 0;
for i = 1:N
    currTargetCalc =  t_sample(:,i)'*(inv(K_avg)*(mu_t_sample - mu_c_sample));
    currClutterCalc = c_sample(:,i)'*(inv(K_avg)*(mu_t_sample - mu_c_sample));
    I_two_rvs_target_response = horzcat(I_two_rvs_target_response, currTargetCalc);
    I_two_rvs_clutter_response = horzcat(I_two_rvs_clutter_response, currClutterCalc);

    if currTargetCalc < eta
        I_two_rvs_prob_miss = I_two_rvs_prob_miss + 1;
    end
    if currClutterCalc > eta
        I_two_rvs_prob_false_alarm = I_two_rvs_prob_false_alarm + 1;
    end
end
I_two_rvs_prob_false_alarm = I_two_rvs_prob_false_alarm / N;
I_two_rvs_prob_miss = I_two_rvs_prob_miss / N;
I_two_rvs_prob_error = I_two_rvs_prob_false_alarm + I_two_rvs_prob_miss;
I_two_rvs_decision_boundary(1,1:N) = eta;

%"dumb" minimization loop to look for MPE(an improvised minimum walking algorithm)
test_eta = [eta];
Prob_errors = [I_two_rvs_prob_error];
for i = 1:5
    degree = 10^-i;
    
    %carry over old candidates
    new_test_eta = test_eta;
    new_Prob_errors = Prob_errors;
    
    for k = 1:length(test_eta)
        currEta = test_eta(k);
        % test back walk
        for j = 1:10
            check_eta = currEta-j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);

            test_back = (sum(I_two_rvs_target_response < check_eta)+sum(I_two_rvs_clutter_response > check_eta))/N;
            new_Prob_errors = horzcat(new_Prob_errors, test_back);
        end

        %test forward walk
        for j = 1:10
            check_eta = currEta+j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);
            test_forward = (sum(I_two_rvs_target_response < check_eta)+sum(I_two_rvs_clutter_response > check_eta))/N; 
            new_Prob_errors = horzcat(new_Prob_errors, test_forward);
        end

        %take our 10 best candidates
        [Prob_errors, errorSort] = sort(new_Prob_errors);
        test_eta = new_test_eta(errorSort);
        if (length(test_eta) > 10)
           Prob_errors = Prob_errors(1:10);
           test_eta = test_eta(1:10);
        end
    end
end
I_two_rvs_min_prob_error = Prob_errors(1);
I_two_rvs_decision_boundary_MPE(1,1:N) = test_eta(1);
sprintf("Probabiliy hunting for two independent r.v.s. found a MPE of %f at a boundary of %f compared to the orginal guess eta of %f's PE of %f",I_two_rvs_min_prob_error, test_eta(1), eta, I_two_rvs_prob_error)


%forming MLR test with tI and cI
K_avg = (KtI+KcI)./2;
mu_t_sample = mu_t(1:3);
mu_c_sample = mu_c(1:3);
eta = 0.5*(mu_t_sample'*inv(K_avg)*mu_t_sample - mu_c_sample'*inv(K_avg)*mu_c_sample);
I_three_rvs_target_response = [];
I_three_rvs_clutter_response = [];
t_sample = [tI{1} tI{2} tI{3}]';
c_sample = [cI{1} cI{2} cI{3}]';
I_three_rvs_prob_false_alarm = 0;
I_three_rvs_prob_miss = 0;
for i = 1:N
    currTargetCalc =  t_sample(:,i)'*(inv(K_avg)*(mu_t_sample - mu_c_sample));
    currClutterCalc = c_sample(:,i)'*(inv(K_avg)*(mu_t_sample - mu_c_sample));
    I_three_rvs_target_response = horzcat(I_three_rvs_target_response, currTargetCalc);
    I_three_rvs_clutter_response = horzcat(I_three_rvs_clutter_response,currClutterCalc);
    
    if currTargetCalc < eta
        I_three_rvs_prob_miss = I_three_rvs_prob_miss + 1;
    end
    if currClutterCalc > eta
        I_three_rvs_prob_false_alarm = I_three_rvs_prob_false_alarm + 1;
    end
end
I_three_rvs_prob_false_alarm = I_three_rvs_prob_false_alarm / N;
I_three_rvs_prob_miss = I_three_rvs_prob_miss / N;
I_three_rvs_prob_error = I_three_rvs_prob_false_alarm + I_three_rvs_prob_miss;
I_three_rvs_decision_boundary(1,1:N) = eta;

%"dumb" minimization loop to look for MPE(an improvised minimum walking algorithm)
test_eta = [eta];
Prob_errors = [I_three_rvs_prob_error];
for i = 1:5
    degree = 10^-i;
    
    %carry over old candidates
    new_test_eta = test_eta;
    new_Prob_errors = Prob_errors;
    
    for k = 1:length(test_eta)
        currEta = test_eta(k);
        % test back walk
        for j = 1:10
            check_eta = currEta-j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);

            test_back = (sum(I_three_rvs_target_response < check_eta)+sum(I_three_rvs_clutter_response > check_eta))/N;
            new_Prob_errors = horzcat(new_Prob_errors, test_back);
        end

        %test forward walk
        for j = 1:10
            check_eta = currEta+j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);
            test_forward = (sum(I_three_rvs_target_response < check_eta)+sum(I_three_rvs_clutter_response > check_eta))/N; 
            new_Prob_errors = horzcat(new_Prob_errors, test_forward);
        end

        %take our 10 best candidates
        [Prob_errors, errorSort] = sort(new_Prob_errors);
        test_eta = new_test_eta(errorSort);
        if (length(test_eta) > 10)
           Prob_errors = Prob_errors(1:10);
           test_eta = test_eta(1:10);
        end
    end
end
I_three_rvs_min_prob_error = Prob_errors(1);
I_three_rvs_decision_boundary_MPE(1,1:N) = test_eta(1);
sprintf("Probabiliy hunting for three independent r.v.s. found a MPE of %f at a boundary of %f compared to the orginal guess eta of %f's PE of %f",I_three_rvs_min_prob_error, test_eta(1), eta, I_three_rvs_prob_error)


figure(19);
h = [];
h(1) = plot(I_one_rvs_target_response, 'DisplayName', 'target reponse');
hold on
h(2) = plot(I_one_rvs_clutter_response, 'DisplayName', 'clutter reponse');
h(3) = plot(I_one_rvs_decision_boundary, 'DisplayName', 'Original Guess Decision Boundary');
h(4) = plot(I_one_rvs_decision_boundary_MPE, 'DisplayName', 'MPE Decision Boundary');
hold off;
title('i.i.d - one r.v MLR discriminator Clutter and Target Responses (using tI and cI)');
axis([0 inf -10 35]);
legend(h);

figure(20);
h = [];
h(1) = plot(I_two_rvs_target_response, 'DisplayName', 'target reponse');
hold on
h(2) = plot(I_two_rvs_clutter_response, 'DisplayName', 'clutter reponse');
h(3) = plot(I_two_rvs_decision_boundary, 'DisplayName', 'Original Guess Decision Boundary');
h(4) = plot(I_two_rvs_decision_boundary_MPE, 'DisplayName', 'MPE Decision Boundary');
hold off;
title('i.i.d - two r.v.s MLR discriminator Clutter and Target Responses (using tI and cI)');
axis([0 inf -10 35]);
legend(h);

figure(21);
h = [];
h(1) = plot(I_three_rvs_target_response, 'DisplayName', 'target reponse');
hold on
h(2) = plot(I_three_rvs_clutter_response, 'DisplayName', 'clutter reponse');
h(3) = plot(I_three_rvs_decision_boundary, 'DisplayName', 'Original Guess Decision Boundary');
h(4) = plot(I_three_rvs_decision_boundary_MPE, 'DisplayName', 'MPE Decision Boundary');
hold off;
title('i.i.d - three r.v.s MLR discriminator Clutter and Target Responses (using tI and cI)');
axis([0 inf -10 35]);
legend(h);


figure(28);
h = [];
h(1) = bar([1 2 3],[I_one_rvs_min_prob_error I_two_rvs_min_prob_error I_three_rvs_min_prob_error], 'DisplayName', 'MPE');
hold on;
h(2) = bar([1 2 3],[I_one_rvs_prob_error I_two_rvs_prob_error I_three_rvs_prob_error], 'DisplayName', 'PE''s under original guess eta', 'FaceColor', 'none', 'EdgeColor', 'red');
hold off;
title('MPE for using 1, 2, and 3 r.v.s for the i.i.d data');
legend(h);
axis([0 inf 0 0.15]);

%C.2.b
%forming MLR test with t{1} and c{1} as they correspond to the largest J
%value
Kc_sample = Kc(1,1);
Kt_sample = Kt(1,1);
mu_t_sample = mu_t_notI(1);
mu_c_sample = mu_c_notI(1);
eta = 0;
one_rvs_target_response = [];
one_rvs_clutter_response = [];
one_rvs_prob_false_alarm = 0;
one_rvs_prob_miss = 0;
for i = 1:N
    currTargetCalc = log(sqrt(det(Kc_sample))) - log(sqrt(det(Kt_sample))) - 0.5*((t{1}(i) - mu_t_sample)'*inv(Kt_sample)*(t{1}(i) - mu_t_sample)) + 0.5*((t{1}(i) - mu_c_sample)'*inv(Kc_sample)*(t{1}(i) - mu_c_sample));
    currClutterCalc = log(sqrt(det(Kc_sample))) - log(sqrt(det(Kt_sample))) - 0.5*((c{1}(i) - mu_t_sample)'*inv(Kt_sample)*(c{1}(i) - mu_t_sample)) + 0.5*((c{1}(i) - mu_c_sample)'*inv(Kc_sample)*(c{1}(i) - mu_c_sample));
    one_rvs_target_response = horzcat(one_rvs_target_response, currTargetCalc);
    one_rvs_clutter_response = horzcat(one_rvs_clutter_response, currClutterCalc);
    
    if currTargetCalc < eta
        one_rvs_prob_miss = one_rvs_prob_miss + 1;
    end
    if currClutterCalc > eta
        one_rvs_prob_false_alarm = one_rvs_prob_false_alarm + 1;
    end
end
one_rvs_prob_false_alarm = one_rvs_prob_false_alarm / N;
one_rvs_prob_miss = one_rvs_prob_miss / N;
one_rvs_prob_error = one_rvs_prob_false_alarm + one_rvs_prob_miss;
one_rvs_decision_boundary(1,1:N) = eta;

%"dumb" minimization loop to look for MPE(an improvised minimum walking algorithm)
test_eta = [eta];
Prob_errors = [one_rvs_prob_error];
for i = 1:5
    degree = 10^-i;
    
    %carry over old candidates
    new_test_eta = test_eta;
    new_Prob_errors = Prob_errors;
    
    for k = 1:length(test_eta)
        currEta = test_eta(k);
        % test back walk
        for j = 1:10
            check_eta = currEta-j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);

            test_back = (sum(one_rvs_target_response < check_eta)+sum(one_rvs_clutter_response > check_eta))/N;
            new_Prob_errors = horzcat(new_Prob_errors, test_back);
        end

        %test forward walk
        for j = 1:10
            check_eta = currEta+j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);
            test_forward = (sum(one_rvs_target_response < check_eta)+sum(one_rvs_clutter_response > check_eta))/N; 
            new_Prob_errors = horzcat(new_Prob_errors, test_forward);
        end

        %take our 10 best candidates
        [Prob_errors, errorSort] = sort(new_Prob_errors);
        test_eta = new_test_eta(errorSort);
        if (length(test_eta) > 10)
           Prob_errors = Prob_errors(1:10);
           test_eta = test_eta(1:10);
        end
    end
end
one_rvs_min_prob_error = Prob_errors(1);
one_rvs_decision_boundary_MPE(1,1:N) = test_eta(1);
sprintf("Probabiliy hunting for one correlated r.v.s. found a MPE of %f at a boundary of %f compared to the orginal guess eta of %f's PE of %f",one_rvs_min_prob_error, test_eta(1), eta, one_rvs_prob_error)


%forming MLR test with t{1}, t{3}, c{1}, and c{3} as they correspond to the two largest J
%values
Kc_sample = Kc([1 3],[1 3]);
Kt_sample = Kt([1 3],[1 3]);
mu_t_sample = mu_t_notI([1 3]);
mu_c_sample = mu_c_notI([1 3]);
eta = 0;
two_rvs_target_response = [];
two_rvs_clutter_response = [];
t_sample = [t{1} t{3}]';
c_sample = [c{1} c{3}]';
two_rvs_prob_false_alarm = 0;
two_rvs_prob_miss = 0;
for i = 1:N
    currTargetCalc = log(sqrt(det(Kc_sample))) - log(sqrt(det(Kt_sample))) - 0.5*((t_sample(:,i) - mu_t_sample)'*inv(Kt_sample)*(t_sample(:,i) - mu_t_sample)) + 0.5*((t_sample(:,i) - mu_c_sample)'*inv(Kc_sample)*(t_sample(:,i) - mu_c_sample));
    currClutterCalc = log(sqrt(det(Kc_sample))) - log(sqrt(det(Kt_sample))) - 0.5*((c_sample(:,i) - mu_t_sample)'*inv(Kt_sample)*(c_sample(:,i) - mu_t_sample)) + 0.5*((c_sample(:,i) - mu_c_sample)'*inv(Kc_sample)*(c_sample(:,i) - mu_c_sample));
    two_rvs_target_response = horzcat(two_rvs_target_response, currTargetCalc);
    two_rvs_clutter_response = horzcat(two_rvs_clutter_response, currClutterCalc);
    
    if currTargetCalc < eta
        two_rvs_prob_miss = two_rvs_prob_miss + 1;
    end
    if currClutterCalc > eta
        two_rvs_prob_false_alarm = two_rvs_prob_false_alarm + 1;
    end
end
two_rvs_prob_false_alarm = two_rvs_prob_false_alarm / N;
two_rvs_prob_miss = two_rvs_prob_miss / N;
two_rvs_prob_error = two_rvs_prob_false_alarm + two_rvs_prob_miss;
two_rvs_decision_boundary(1,1:N) = eta;

%"dumb" minimization loop to look for MPE(an improvised minimum walking algorithm)
test_eta = [eta];
Prob_errors = [two_rvs_prob_error];
for i = 1:5
    degree = 10^-i;
    
    %carry over old candidates
    new_test_eta = test_eta;
    new_Prob_errors = Prob_errors;
    
    for k = 1:length(test_eta)
        currEta = test_eta(k);
        % test back walk
        for j = 1:10
            check_eta = currEta-j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);

            test_back = (sum(two_rvs_target_response < check_eta)+sum(two_rvs_clutter_response > check_eta))/N;
            new_Prob_errors = horzcat(new_Prob_errors, test_back);
        end

        %test forward walk
        for j = 1:10
            check_eta = currEta+j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);
            test_forward = (sum(two_rvs_target_response < check_eta)+sum(two_rvs_clutter_response > check_eta))/N; 
            new_Prob_errors = horzcat(new_Prob_errors, test_forward);
        end

        %take our 10 best candidates
        [Prob_errors, errorSort] = sort(new_Prob_errors);
        test_eta = new_test_eta(errorSort);
        if (length(test_eta) > 10)
           Prob_errors = Prob_errors(1:10);
           test_eta = test_eta(1:10);
        end
    end
end
two_rvs_min_prob_error = Prob_errors(1);
two_rvs_decision_boundary_MPE(1,1:N) = test_eta(1);
sprintf("Probabiliy hunting for two correlated r.v.s. found a MPE of %f at a boundary of %f compared to the orginal guess eta of %f's PE of %f",two_rvs_min_prob_error, test_eta(1), eta, two_rvs_prob_error)


%forming MLR test with t and c
Kc_sample = Kc;
Kt_sample = Kt;
mu_t_sample = mu_t_notI;
mu_c_sample = mu_c_notI;
eta = 0;
three_rvs_target_response = [];
three_rvs_clutter_response = [];
t_sample = [t{1} t{2} t{3}]';
c_sample = [c{1} c{2} c{3}]';
three_rvs_prob_false_alarm = 0;
three_rvs_prob_miss = 0;
for i = 1:N
    currTargetCalc = log(sqrt(det(Kc_sample))) - log(sqrt(det(Kt_sample))) - 0.5*((t_sample(:,i) - mu_t_sample)'*inv(Kt_sample)*(t_sample(:,i) - mu_t_sample)) + 0.5*((t_sample(:,i) - mu_c_sample)'*inv(Kc_sample)*(t_sample(:,i) - mu_c_sample));
    currClutterCalc = log(sqrt(det(Kc_sample))) - log(sqrt(det(Kt_sample))) - 0.5*((c_sample(:,i) - mu_t_sample)'*inv(Kt_sample)*(c_sample(:,i) - mu_t_sample)) + 0.5*((c_sample(:,i) - mu_c_sample)'*inv(Kc_sample)*(c_sample(:,i) - mu_c_sample));
    three_rvs_target_response = horzcat(three_rvs_target_response, currTargetCalc);
    three_rvs_clutter_response = horzcat(three_rvs_clutter_response, currClutterCalc);

    if currTargetCalc < eta
        three_rvs_prob_miss = three_rvs_prob_miss + 1;
    end
    if currClutterCalc > eta
        three_rvs_prob_false_alarm = three_rvs_prob_false_alarm + 1;
    end
end
three_rvs_prob_false_alarm = three_rvs_prob_false_alarm / N;
three_rvs_prob_miss = three_rvs_prob_miss / N;
three_rvs_prob_error = three_rvs_prob_false_alarm + three_rvs_prob_miss;
three_rvs_decision_boundary(1,1:N) = eta;

%"dumb" minimization loop to look for MPE(an improvised minimum walking algorithm)
test_eta = [eta];
Prob_errors = [three_rvs_prob_error];
for i = 1:5
    degree = 10^-i;
    
    %carry over old candidates
    new_test_eta = test_eta;
    new_Prob_errors = Prob_errors;
    
    for k = 1:length(test_eta)
        currEta = test_eta(k);
        % test back walk
        for j = 1:10
            check_eta = currEta-j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);

            test_back = (sum(three_rvs_target_response < check_eta)+sum(three_rvs_clutter_response > check_eta))/N;
            new_Prob_errors = horzcat(new_Prob_errors, test_back);
        end

        %test forward walk
        for j = 1:10
            check_eta = currEta+j*degree;
            new_test_eta = horzcat(new_test_eta, check_eta);
            test_forward = (sum(three_rvs_target_response < check_eta)+sum(three_rvs_clutter_response > check_eta))/N; 
            new_Prob_errors = horzcat(new_Prob_errors, test_forward);
        end

        %take our 10 best candidates
        [Prob_errors, errorSort] = sort(new_Prob_errors);
        test_eta = new_test_eta(errorSort);
        if (length(test_eta) > 10)
           Prob_errors = Prob_errors(1:10);
           test_eta = test_eta(1:10);
        end
    end
end
three_rvs_min_prob_error = Prob_errors(1);
three_rvs_decision_boundary_MPE(1,1:N) = test_eta(1);
sprintf("Probabiliy hunting for three correlated r.v.s. found a MPE of %f at a boundary of %f compared to the orginal guess eta of %f's PE of %f",three_rvs_min_prob_error, test_eta(1),eta, three_rvs_prob_error)



figure(22);
h = [];
h(1) = plot(one_rvs_target_response, 'DisplayName', 'target reponse');
hold on
h(2) = plot(one_rvs_clutter_response, 'DisplayName', 'clutter reponse');
h(3) = plot(one_rvs_decision_boundary, 'DisplayName', 'Original Guess Decision Boundary');
h(4) = plot(one_rvs_decision_boundary_MPE, 'DisplayName', 'MPE Decision Boundary');
hold off;
title('correlated - one r.v MLR discriminator Clutter and Target Responses (using t and c)');
axis([0 inf -40 40]);
legend(h);



figure(23);
h = [];
h(1) = plot(two_rvs_target_response, 'DisplayName', 'target reponse');
hold on
h(2) = plot(two_rvs_clutter_response, 'DisplayName', 'clutter reponse');
h(3) = plot(two_rvs_decision_boundary, 'DisplayName', 'Original Guess Decision Boundary');
h(4) = plot(two_rvs_decision_boundary_MPE, 'DisplayName', 'MPE Decision Boundary');
hold off;
title('correlated - two r.v.s MLR discriminator Clutter and Target Responses (using t and c)');
axis([0 inf -40 40]);
legend(h);

figure(24);
h = [];
h(1) = plot(three_rvs_target_response, 'DisplayName', 'target reponse');
hold on
h(2) = plot(three_rvs_clutter_response, 'DisplayName', 'clutter reponse');
h(3) = plot(three_rvs_decision_boundary, 'DisplayName', 'Original Guess Decision Boundary');
h(4) = plot(three_rvs_decision_boundary_MPE, 'DisplayName', 'MPE Decision Boundary');
hold off;
title('correlated - three r.v.s MLR discriminator Clutter and Target Responses (using t and c)');
axis([0 inf -40 40]);
legend(h);

figure(29);
h = [];
h(1) = bar([1 2 3],[one_rvs_min_prob_error two_rvs_min_prob_error three_rvs_min_prob_error], 'DisplayName', 'MPE');
hold on;
h(2) = bar([1 2 3],[one_rvs_prob_error two_rvs_prob_error three_rvs_prob_error], 'DisplayName', 'PE''s under original guess eta', 'FaceColor', 'none', 'EdgeColor', 'red');
hold off;
title('MPE for using 1, 2, and 3 r.v.s for the correlated data');
legend(h);
axis([0 inf 0 0.15]);