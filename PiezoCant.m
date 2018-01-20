%% ------Energy Harvesting using Piezo Material with Cantilever Structures------

%  This promgram intend to use the Continuous Vibration Equation for estimate 
%  the power generation of a cantilever fully covered by piezoelectric material
%  The Math background is provided by the paper "A Distributed Parameter 
%  Electromechanical Model For Cantilevered Piezoelectric Energy Harvesters" by
%  A.Erturk and D.J.Inman.


clear ; close all; clc


%----------------------------INPUT PARAMETERS----------------------------------%

% This simulation can range diffente input parameters and then plot it. Choose
% one parameter and replace by MRANGE(c,2) variable. The MRANGE(c,2) will change
% According to the RANGE Variable bellow

% RANGE MODULE

RANGE = 0.1:0.1:0.5;         % Range for used for change the parameters (MAX: 8)
VRANGE = 1:length(RANGE);    % Vector with same size of the RANGE;
MRANGE = [VRANGE',RANGE'];   % Matrix with The the Index of RANGE

F = 1:0.1:100;               % Frequency Range (Hz)

RANGEV0 = ones(length(F),length(RANGE)); % Matrix for Voltage simulation data
RANGEI0 = ones(length(F),length(RANGE)); % Matrix for Current simulation data
RANGEP0 = ones(length(F),length(RANGE)); % Matrix for Power simulation data


% MECHANICAL PARAMETER INPUT

% Beam Support

for c = VRANGE         
 
 leB = MRANGE(c,2);    % Lenght of the Beam (m)
 wiB = 0.02;           % Widht of the Beam (m)
 thB = 0.0005;         % Thickness of the Beam (m)
 yoB = 100000000000;   % Young's Modulos of the Beam kgf/m2
 deB = 7165;           % Material Density (kg/m3)


% Piezoelectric material

 leP = leB;            % Lenght of the Piezo (m) OBS: SHOULD BE SAME AS THE BEAM
 wiP = wiB;            % Widht of the Piezo (m)
 thP = 0.0004;         % Thickness of the Piezo (m)
 deP = 7800;           % Material Density of Piezo (kg/m3)
 yoP = 66000000000;    % Young's Modulos of the Piezo


% ELECTRICAL PARAMETER INPUT

% Piezoelectric material

 d31 = -190*(10^-12);  % Piezoelectric Constant at 31 Direction 
 e33 = 15.93*(10^-9);  % Piezoelectric Dielectric Constant (F/m)



% BASE HARMONIC VIBRATION

 Y0 = 0.000001;        % Peak Displacement of the base (m)
 f = F;                % Frequency Range (Hz)
 w = 2*pi()*f;         % Angular frequency (Hz)

% ---------------------------FIRST CALCULATIONS--------------------------------%




% CALCULATING THE MASS DENSINTY

% The mass density is the total mass of the system dived by the total length of
% the Beam

 volB = thB*wiB*leB;
 volP = thP*wiP*leP;
 maB = deB*volB;
 maP = deP*volP;

 massden = (maB+maP)/leB; % (kg/m)

% LOAD RESISTOR

% Resistor conected to the Piezo

 rl = 1000000 ;% Ohms


% POSITION RELATIVE TO NEUTRAL AXIS

% The procedure of finding the position of the neutral axis of a composite cross
% section is transform the cross section to a homogeneous cross section of 
% single Young’s modulus (Described in Timoshenko and Young Reference)

 nY = yoB/yoP;                                                                   % Ratio of Young's Modulos

 hpa = ((thP ^ 2) + (2*nY*thP*thB) + (nY*(thB ^ 2)))/(2*(thP + (nY*thB)));       % Distance from the top of the Piezo Layer to the Neutral Axis

 hsa = ((thP ^ 2) + (2*thP*thB) + (nY*(thB ^ 2)))/(2*(thP + (nY*thB)));          % Distance from the Bottom of the Piezo Layer to the Neutral Axis

 hpc = (nY*thB*(thP + thB))/(2*(thP + (nY*thB)));                                % Distance from the Center of the Piezo Layer to the Neutral Axis


 ha = -hsa;                                                                      % Position of the bottom of the substructure layer from the neutral axis
 hb = hpa - thP;                                                                 % Position of the bottom of the Piezo layer from the neutral axis
 hc = hpa;                                                                       % Position of the top of the  layer from the neutral axis


% BENDING STIFFNESS 

 yi = wiB*(((yoB*((hb^3)-(ha^3))) + (yoP*((hc^3) - (hb^3))))/3);


% COUPLING TERM

 coup = - (((yoP*d31*wiB) / ((2*thP))*((hc^2) -(hb^2))));


% UNDAMPED NATURAL FREQUENCY OF THE BEAM

% The sigma constant for Cantilever Structuree is provide by the literature.
% The Solution for sigma is solving the equation 1 + cos (sig) cosh (sig) = 0
% For the three first modes we have: 

 sig1 = 1.875104069;
 sig2 = 4.694091133;
 sig3 = 7.854757438;

 sig = [sig1;sig2;sig3];

 wr1 = (sig1^2) * ((yi / (massden*(leB^4)))^0.5);
 wr2 = (sig2^2) * ((yi / (massden*(leB^4)))^0.5);
 wr3 = (sig3^2) * ((yi / (massden*(leB^4)))^0.5);

 wr = [wr1;wr2;wr3];

% EIGENFUNCTIONS

 ho = ((cosh(sig) + cos(sig))./(sin(sig) + sinh(sig)));
 x = leB;
 var = ((sig*x)./leB);
 eigeinf = ((1/(massden*leB))^0.5).*(cosh(var) - cos(var) - (ho.*(sinh(var)-sin(var))));

 
% BENDING SLOPE (DERIVATE OF EIGENFUNCTIONS)

 bendS = ((1/(massden*leB))^0.5)*(((sig./leB).* sinh(var)) + ((sig./leB).* sin(var)) - (ho.*(((sig./leB).*cosh(var))-((sig./leB).* cos(var)))));


% INTEGRATION OF THE EIGENFUNCTIONS

 inteia = ((1/(massden*leB))^0.5)*((sinh(var)./(sig./leB)) - (sin(var)./(sig./leB)) - (ho.*((cosh(var)./(sig./leB))+(cos(var)./(sig./leB)))));
 intei0 = ((1/(massden*leB))^0.5)*((sinh(0)./(sig./leB)) - (sin(0)./(sig./leB)) - (ho.*((cosh(0)./(sig./leB))+(cos(0)./(sig./leB)))));
 intei = inteia - intei0;


% TIME CONSTANT OF THE PIEZO REPRESENTATION CIRCUIT

 tc = (rl*e33*wiP*leP) / thP;


% MODAL CONSTANT TERM

 mco = -(((d31*yoP*hpc*thP)/(e33*leB)) * bendS);


% MODAL COUPLING TERM

 mc = coup * bendS;


% MODAL DAMPING (ASSUMING DAMPING RATIOS ZETA 1 = 0.01 AND ZETA 2 = 0.013)

 ca = 0.654961133;
 csI = 1.216501665*(10^-6);
 zeta = ((csI*wr)/(2*yi))+ (ca./(2*wr*massden));


% VOLTAGE ESTIMATION

 modes = 3;
 fs = size(f);

 numV = ones(modes,fs(2));
 denV = ones(modes,fs(2));

 for b = 1:modes
 
 numV(b,:) = ((-i*massden*w*mco(b)*intei(b))./(((wr(b)^2)-(w.^2)) + (i*2*zeta(b)*wr(b)*w))); 
 denV(b,:) = ((i*w*mc(b)*mco(b))./(((wr(b)^2)-(w.^2)) + (i*2*zeta(b)*wr(b)*w)));
 
 endfor

 num = numV(1,:) + numV(2,:) + numV(3,:);
 den = denV(1,:) + denV(2,:) + denV(3,:);

 denV2 =((1 + (i*w*tc))/tc);

 v = num./(den+denV2);
 V0 = abs(v);
 V = V0.*((w.^2)*Y0);
 RANGEV0(:,MRANGE(c,1)) = V0';
 RANGEV(:,MRANGE(c,1)) = V';


% CURRENT ESTIMATION

 j = (num./((den+denV2)*rl));
 I0 = abs(j);
 I = I0.*((w.^2)*Y0);
 RANGEI0(:,MRANGE(c,1)) = I0';
 RANGEI(:,MRANGE(c,1)) = I';

% POWER ESTIMATION

 p = v.*j;
 P0 = abs(p);
 P = P0.*(((w.^2)*Y0).^2);
 RANGEP0(:,MRANGE(c,1)) = P0';
 RANGEP(:,MRANGE(c,1)) = P';
 
endfor


% PLOTTING DATA

% Input the legend, title, and lables for Axis. Assign the Range variable
% VOLTAGE, CURRENT or POWER for plot (Use RANGEV0, RANGEI0 or RANGEP0 respectively). 

xlab = "Frequency (Hz)";
ylab = "FRF Voltage (V / Base Acceleration) ";
tit =  "FRF Voltage Variation with the Length(m)";
legS = "Length ";
Range = RANGEV0;

figure(1);
 labels = {};
 colororder = get (gca, "colororder");
 for c = VRANGE
   h = plot (f,Range(:,c));
   hold on;
   set (h, "color", colororder(c,:));
   labels = {labels{:}, [legS, num2str(MRANGE(c,2))]};
 endfor
 hold off;
 hx = xlabel({xlab},"fontsize",12);

 ylabel(ylab,"fontsize", 12);
 grid on;
 title (tit,"fontsize",12);
 legend (labels);
 print -djpg figure1.jpg;
 

% EXPORTING DATA 

% You can export data in thge .mat format. Just remove the comments bellow

% FData = [f',RANGEV0,RANGEI0,RANGEP0];
% save CantileverFRF.mat FData;
% FData = [f',RANGEV,RANGEI,RANGEP];
% save CantileverRealV.mat FData;