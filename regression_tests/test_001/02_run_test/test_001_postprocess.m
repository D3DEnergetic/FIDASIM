% Postprocess test_001 results:

clear all
close all
clc

% Flags:
save_figure = 0;
save_data = 0;

% Define directories:
output_dir = "./";
input_dir = "./";

%% GET data:

test_file = fullfile(input_dir,"test_001.h5");
info = h5info(test_file);

% Grids:
x = h5read(test_file,'/x');
dene = h5read(test_file,'/dene');
T = h5read(test_file,'/temperature'); % [eV]
eb = h5read(test_file,'/eb');

% Density profiles:
denn = h5read(test_file,'/denn');
denn_n = h5read(test_file,'/denn_n');

% Dimensions:
nden = numel(dene);
neb = numel(eb);
nt = numel(T);

%% COMPUTE mean free path:
for nn = 1:nden
    for ee = 1:neb
        for tt = 1:nt
            [Ls,y] = compute_mfp(x,denn,nn,ee,tt);
            mfp(nn,ee,tt) = Ls;
            yarr{nn,ee,tt} = y;
        end
    end
end

%% PLOT data:

font_size.ax = 14;
font_size.title = 16;
font_size.label = 16;
font_size.legend = 14;

line_color = {"k","r","g","bl","m","c"};

% Loop over density vales:
for nn = 1:nden
    figure('color','w')
    hold on
    box on
    ymax = 0;
    set(gca,'FontSize',font_size.ax)

    % Loop over energy:
    for ee = 1:neb
        yval = (squeeze(mfp(nn,ee,:)))';
        ymax = max(max(yval),ymax);

        % Plot vs temperature:
        hls(ee) = plot(T*1e-3,yval,line_color{ee},'linewidth',3);
        leg_str{ee} = "$n_e$: " + sprintf('%.1e',dene(nn)) + " [cm$^{-3}$], " + ...
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
        figure_name = "FIDASIM" + "_MFP_" + nn;
        save_figure_to_file("./figures",figure_name)
    end

end

%% APPEND data to test file:

if save_data
    % Create new dataset:
    try
        h5create(test_file,"/mfp",size(mfp),"Datatype","double")
    end

    % Write data:
    h5write(test_file,"/mfp",mfp);
    
    % Add attributes:
    h5writeatt(test_file,"/mfp", "description", "Mean free path computed from log-linear fit of total neutral density from COLRAD");
    h5writeatt(test_file,"/mfp", "units", "cm");
    h5writeatt(test_file,"/mfp", "dimensions", "[dene, eb, temperature]");
end

%% PLOT 1D profile example:

% Select data:
nn = 1;
ee = 1;
tt = 1;

% Extract values
ne  = dene(nn);   % [cm^-3]
TkeV = T(tt)/1e3; % [keV]
EbkeV = eb(ee)/1e3; % [keV/amu]

% Format density in scientific notation (mantissa × 10^{exp})
exp_ne = floor(log10(ne));
mant_ne = ne / 10^exp_ne;

title_str = sprintf([ ...
  '$n_e = %.2f\\times 10^{%d}\\,\\mathrm{cm}^{-3},\\ ' ...
  'T = %.3g\\,\\mathrm{keV},\\ ' ...
  'E_b = %.3g\\,\\mathrm{keV/amu}$'], ...
  mant_ne, exp_ne, TkeV, EbkeV);

% Normalize profile:
y = denn(:,nn,ee,tt);
y = y/max(y);

figure('color','w')
hold on
box on
grid on
set(gca,'FontSize',font_size.ax)
plot(x(:,nn),log10(y),'LineWidth',3)
set(gca,"yscale",'lin')
title(title_str,"interpreter","latex","FontSize",font_size.title)
xlabel("x [cm]","interpreter","latex","FontSize",font_size.label)
ylabel("log$_{10}(n_n/n_{n0})$","interpreter","latex","FontSize",1.5*font_size.label)

if save_figure
    figure_name = "Example_neutral_density_profile";
    save_figure_to_file("./figures",figure_name)
end

%% Function: compute_mfp
function [Ls,y] = compute_mfp(x,denn,nn,ee,tt)

dx = x(2,nn) - x(1,nn);
z = denn(:,nn,ee,tt)/max(denn(:,nn,ee,tt));
dlogz = diff(log(z))/dx;
y = -1./dlogz;
Ls = mean(y);

end