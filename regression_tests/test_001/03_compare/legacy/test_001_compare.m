% Compare reference and test results
clear all
close all
clc

% Flags:
save_figure = 1;
save_data = 0;

% Acceptance criterion:
mfp_relative_error_tolerance = 0.12; % Maximum allowed relative error (12%).

%% GET reference data:
source_dir = "../01_reference/";
file_name = "reference.h5";
file_path = fullfile(source_dir,file_name);
info{1} = h5info(file_path);

% Read data in
data{1}.T   = h5read(file_path,"/T");
data{1}.eb  = h5read(file_path,"/eb");
data{1}.ne  = h5read(file_path,"/ne");
data{1}.mfp = h5read(file_path,"/mfp");
data{1}.ab  = h5read(file_path,"/ab");

%% GET test data:
source_dir = "../02_run_test/";
file_name = "test_001.h5";
file_path = fullfile(source_dir,file_name);
info{2} = h5info(file_path);

% Read data in
data{2}.T   = h5read(file_path,"/temperature");
data{2}.eb  = h5read(file_path,"/eb");
data{2}.ne  = h5read(file_path,"/dene");
data{2}.mfp = h5read(file_path,"/mfp");
data{2}.ab  = h5read(file_path,"/beam_mass");
data{2}.x   = h5read(file_path,"/x");
data{2}.denn = h5read(file_path,"/denn");

font_size.ax = 14;
font_size.title = 16;
font_size.label = 16;
font_size.legend = 14;

dene = data{1}.ne;
eb = data{1}.eb;
nden = numel(data{1}.ne);
neb = numel(data{1}.eb);
nt = numel(data{1}.T);

%% PLOT data comparison:

line_color = {"k","r","g","bl","m","c"};

% Loop over density vales:
for nn = 1:nden
    figure('color','w')
    set(gcf,"position",[675, 180, 570, 750]) % [left bottom width height]
    hold on
    box on
    ymax = 0;
    set(gca,'FontSize',font_size.ax)

    % Loop over energy:
    for ee = 1:neb

        % ref:
        yval = (squeeze(data{1}.mfp(nn,ee,:)))';
        ymax = max(max(yval),ymax);
        hls(ee) = plot(data{1}.T*1e-3,yval,line_color{ee},'linewidth',3);

        % test:
        yval = (squeeze(data{2}.mfp(nn,ee,:)))';
        ymax = max(max(yval),ymax);
        plot(data{2}.T*1e-3,yval,line_color{ee},'linewidth',2,'LineStyle','--');


        leg_str{ee} = "(REF) $n_e$: " + sprintf('%.1e',dene(nn)) + " [cm$^{-3}$], " + ...
                      "$E_b$: " + sprintf('%.1e',eb(ee)*1e-3) + " [keV]";
    end
    grid on
    
    hleg = legend(hls,leg_str);
    set(hleg,'Interpreter','latex','FontSize',font_size.legend)
    set(hleg,'Location','eastoutside')
    set(hleg,'Location','southoutside')
    
    title("Mean free path $\lambda_s$ [cm]",'Interpreter','latex',FontSize=font_size.title)
    ylabel("$\lambda_s$ [cm]", 'Interpreter','latex','FontSize',font_size.label)
    ylim([0,1.2]*ymax)
    xlabel("T [keV]", 'Interpreter','latex','FontSize',font_size.label)          
    
    if save_figure
        figure_name = "comparison_" + "_MFP_" + nn;
        save_figure_to_file("./figures",figure_name)
    end

end

%% COMPUTE relative error:
mfp_test = data{2}.mfp;
mfp_ref = data{1}.mfp;

rel_err = abs(mfp_test - mfp_ref) ./ abs(mfp_ref);

max_err  = max(rel_err(:));
mean_err = mean(rel_err(:),'omitnan');
std_err = std(rel_err(:),'omitnan');

regression_passed = isfinite(max_err) && ...
                    max_err <= mfp_relative_error_tolerance;
if regression_passed
    regression_status = 'PASS';
else
    regression_status = 'FAIL';
end

%% PLOT relative error:

report = sprintf([ ...
    '==============================================\n' ...
    ' OVERALL REGRESSION STATUS: %s\n' ...
    ' Maximum MFP relative error [%%]: %.2f\n' ...
    ' Acceptance tolerance [%%]       : %.2f\n' ...
    ' Acceptance check: %.2f%% <= %.2f%%\n' ...
    '==============================================\n' ...
    '\n' ...
    '----------------------------------------------\n' ...
    ' Mean Free Path Regression Test Summary\n' ...
    '----------------------------------------------\n' ...
    '  Max relative error [%%]   : %.2f\n' ...
    '  Mean relative error [%%]  : %.2f\n' ...
    '  Std. relative error [%%]  : %.2f\n' ...
    '----------------------------------------------\n' ...
    '\n'], ...
    regression_status, max_err*1e2, mfp_relative_error_tolerance*1e2, ...
    max_err*1e2, mfp_relative_error_tolerance*1e2, ...
    max_err*1e2, mean_err*1e2, std_err*1e2);

fprintf('%s', report);

if save_data
    outfile = 'mfp_regression_report.txt';
    fid = fopen(outfile, 'w');
    if fid < 0
        error('Could not open %s for writing', outfile);
    end
    fprintf(fid, '%s', report);
    fclose(fid);
end

% Loop over density vales:
for nn = 1:nden
    figure('color','w')
    hold on
    box on
    ymax = 0;
    set(gca,'FontSize',font_size.ax)

    % Loop over energy:
    for ee = 1:neb

        yval = (squeeze(rel_err(nn,ee,:)*1e2))';
        ymax = max(max(yval),ymax);
        hls(ee) = plot(data{1}.T*1e-3,yval,line_color{ee},'linewidth',3);

        leg_str{ee} = "$n_e$: " + sprintf('%.1e',dene(nn)) + " [cm$^{-3}$], " + ...
                      "$E_b$: " + sprintf('%.1e',eb(ee)*1e-3) + " [keV]";
    end
    grid on
    
    hleg = legend(hls,leg_str);
    set(hleg,'Interpreter','latex','FontSize',font_size.legend)
    set(hleg,'Location','eastoutside')
    set(hleg,'Location','southoutside')
    
    title("Relative error [$\%$]",'Interpreter','latex',FontSize=font_size.title)
    ylabel("rel. err. [$\%$]", 'Interpreter','latex','FontSize',font_size.label)
    ylim([0,mfp_relative_error_tolerance]*1e2)
    xlabel("T [keV]", 'Interpreter','latex','FontSize',font_size.label)          
    
    if save_figure
        figure_name = "comparison_" + "_relative_error_" + nn;
        save_figure_to_file("./figures",figure_name)
    end

end

%% PLOT profile with largest relative error:

% Options:
y_scale = 'lin'; % lin or log

% Get indices of largest error:
[~, idx] = max(rel_err(:), [], 'omitnan');
[nn,ee,tt] = ind2sub(size(rel_err), idx);

% nn = 2;
% ee = 2;
% tt = 10;

% Extract values
ne  = data{1}.ne(nn);   % [cm^-3]
TkeV = data{1}.T(tt)/1e3; % [keV]
EbkeV = data{1}.eb(ee)/1e3; % [keV/amu]

% Format density in scientific notation (mantissa × 10^{exp})
exp_ne = floor(log10(ne));
mant_ne = ne / 10^exp_ne;

title_str = sprintf([ ...
  '$n_e = %.2f\\times 10^{%d}\\,\\mathrm{cm}^{-3},\\ ' ...
  'T = %.3g\\,\\mathrm{keV},\\ ' ...
  'E_b = %.3g\\,\\mathrm{keV/amu}$'], ...
  mant_ne, exp_ne, TkeV, EbkeV);

% Get profiles:
x = data{2}.x(:,nn);
y = data{2}.denn(:,nn,ee,tt);
denn_test = y/max(y);
lambda_s = data{1}.mfp(nn,ee,tt);
denn_ref = exp(-x/lambda_s);

% Plot data comparison:
figure('color','w')
box on
hold on
set(gca,'FontSize',font_size.ax)
hp(1) = plot(x,denn_ref,'k','LineWidth',3);
legend_str{1} ="REF";
hp(2) = plot(x,denn_test,'k--','LineWidth',3);
legend_str{2} ="TEST";
xlim([0,3*lambda_s])
set(gca,'YScale',y_scale)

hleg = legend(hp,legend_str);
set(hleg,"FontSize",font_size.legend,'Interpreter','latex')

xlabel("x [cm]","interpreter","latex","FontSize",font_size.label)
switch y_scale
    case "lin"
        ylabel("$n_n/n_{n0}$","interpreter","latex","FontSize",1.5*font_size.label)
    case "log"
        ylabel("log$_{10}(n_n/n_{n0})$","interpreter","latex","FontSize",1.5*font_size.label)
end
title(title_str,"interpreter","latex","FontSize",font_size.title)


if save_figure
    figure_name = "max" + "_relative_error_profiles";
    save_figure_to_file("./figures",figure_name)
end

if ~regression_passed
    error('test_001:RegressionFailed', ...
        ['Test 001 failed: maximum MFP relative error %.2f%% exceeds ' ...
         'the %.2f%% tolerance.'], ...
        max_err*1e2, mfp_relative_error_tolerance*1e2);
end

return

% Look at the power deposition strength:
figure
hold on
plot(x(1:end-1),-diff(denn_ref)./diff(x),'k','LineWidth',3);
plot(x(1:end-1),-diff(denn_test)./diff(x),'k--','LineWidth',3);
