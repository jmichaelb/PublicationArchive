function [ECStruct,Cout,ea]=mkStrGlaucophane(Cflg)
% L. Bezacier, B. Reynard, J.D. Bass, J. Wang, D. Mainprice
% Elasticity of glaucophane, seismic velocities and anisotropy of the subducted
%  oceanic crust  Tectonophysics 494 (2010) 201?210

ECStruct.Data.name='Glaucophane from Bezacier et al (2010)';

% the default is to use the published direction cosines for each velocity
% measurment.  Changing this to 0 allows euler angles to be used
ECStruct.Data.dcosflg=1;

% the data were taken on three separate slices of the crystal
ECStruct.Data.nsamp=3;

% set symmetry and density
ECStruct.Data.sym='m';
ECStruct.Data.rho=3.07;

% starting estimates for the moduli (cyclic order)
% use published moduli or return default values
if strcmp(Cflg,'p')
    ECStruct.Data.Cguess=[122.3
    45.7
    37.2
    2.3
    231.5
    74.9
    -4.8
    254.6
    -23.7
    79.6
    8.9
    52.8
    51.2];
else
ECStruct.Data.Cguess=[100    
                       50  
                       50  
                        0   
                      100   
                       50 
                        0  
                      100   
                       50  
                       50  
                        0  
                       50   
                       50   ];
end
                     


% set "trust region"
ECStruct.Data.Trust.constants=[50   300
                       10   100
                       10   100
                      -10   100
                       50   400
                       10   100
                      -30    30
                       50   400
                      -30    50
                       15   100
                      -50    50
                       10   100
                       10   100
];
ECStruct.Data.Trust.eulerangles=[20;20;20]*ones(1,ECStruct.Data.nsamp);

% experimental data for velocities copied from paper
velocities=  [ NaN    4.4170    7.1140
    4.2270       NaN    7.7740
    4.1350    4.8370    8.3360
    4.0450    5.0580    8.6340
    4.0670    5.2080    8.6190
    4.1450    5.1980    8.2580
    4.2720    5.0770    7.7060
    4.3340    4.8140    7.0510
       NaN    4.3650    6.6490
       NaN    4.2700    6.4630
       NaN    4.2300    6.3640
       NaN    4.3940    6.5900
       NaN    4.5010    7.1840
       NaN    4.2890    6.2340
    4.1780    4.5530    6.4800
    4.6950       NaN    6.8050
    4.6970       NaN    7.6020
    4.4040       NaN    8.5280
    4.0920       NaN    8.9900
    3.9390       NaN    9.1490
    4.0320       NaN    8.9120
    4.2730    4.6780    8.3500
    4.2660    4.7990    7.4040
    4.1590    4.8660    6.5930
       NaN    4.2960    6.3330
    3.9950    5.1220    8.5490
    4.0340    5.0760    8.6420
    4.1010    5.1180    8.7160
    4.1260    5.2210    8.7910
    4.0070    5.2240    8.9780
       NaN    5.0520    9.1360
    3.8960    4.9230    9.2300
    3.8710    4.8900    9.1430
    3.9740    5.0270    8.8900
    4.0340    5.1770    8.6330
    4.1430    5.2440    8.4640
    4.0690    5.2000    8.4350
    4.0240    5.1140    8.5410];

%direction cosines copied from paper
dcos=[0.8070    0.5270    0.2650
    0.6300    0.7120    0.3090
    0.4110    0.8500    0.3310
    0.1630    0.9290    0.3320
   -0.0960    0.9460    0.3110
   -0.3480    0.8980    0.2690
   -0.5780    0.7890    0.2110
   -0.7680    0.6260    0.1370
   -0.9060    0.4210    0.0530
   -0.9820    0.1870   -0.0350
   -0.9910   -0.0590   -0.1240
   -0.9310   -0.3010   -0.2050
   -0.8080   -0.5220   -0.2740
    0.9920   -0.1140   -0.0510
    0.8980    0.0310    0.4390
    0.7560    0.0980    0.6470
    0.5630    0.1590    0.8110
    0.3330    0.2090    0.9190
    0.0810    0.2470    0.9660
   -0.1770    0.2710    0.9460
   -0.4240    0.2790    0.8620
   -0.6420    0.2690    0.7180
   -0.8180    0.2400    0.5230
   -0.9370    0.1930    0.2920
   -0.9910    0.1310    0.0420
    0.2860    0.9570   -0.0390
    0.2190    0.9510    0.2170
    0.1350    0.8790    0.4570
    0.0410    0.7460    0.6650
   -0.0570    0.5600    0.8260
   -0.1490    0.3360    0.9300
   -0.2300    0.0890    0.9690
   -0.2930   -0.1640    0.9420
   -0.3350   -0.4050    0.8510
   -0.3530   -0.6180    0.7020
   -0.3480   -0.7890    0.5070
   -0.3210   -0.9060    0.2780
   -0.2730   -0.9610    0.0300];

% euler orientations of the crystal in lab coordinates determined by taking 
% cross product of direction cosines

ECStruct.Data.eulerangles=[
  253.16  264.22  343.93
   19.75   75.09   75.77
  141.68   86.97   87.69];

fac=.008;  % uncertainty estimate based on published misfit of 44 m/s

%Load the structure
ECStruct.Data.sample(1).name='Glaucophane 1';
ECStruct.Data.sample(1).BWvelocities= velocities(1:13,:);
ECStruct.Data.sample(1).BWangles= [
     0
    15
    30
    45
    60
    75
    90
   105
   120
   135
   150
   165
   179];
ECStruct.Data.sample(1).BWuncertainties=fac*ECStruct.Data.sample(1).BWvelocities;
ECStruct.Data.sample(1).dcos= dcos(1:13,:);
 
ECStruct.Data.sample(2).name='Glaucophane 2';     
ECStruct.Data.sample(2).BWvelocities= velocities(14:25,:);
ECStruct.Data.sample(2).dcos= dcos(14:25,:);
ECStruct.Data.sample(2).BWangles= [
     0
    30
    45
    60
    75
    90
   105
   120
   135
   150
   165
   179];
ECStruct.Data.sample(2).BWuncertainties=fac*ECStruct.Data.sample(2).BWvelocities;
  
ECStruct.Data.sample(3).name='Glaucophane 3';
ECStruct.Data.sample(3).BWvelocities= velocities(26:38,:);
ECStruct.Data.sample(3).dcos= dcos(26:38,:);
ECStruct.Data.sample(3).BWangles=[
     0
    15
    30
    46
    61
    76
    91
   106
   121
   136
   151
   166
   179];
ECStruct.Data.sample(3).BWuncertainties=fac*ECStruct.Data.sample(3).BWvelocities;

%   Set the opts variables
ECStruct.opts.eulerangles=ECStruct.Data.eulerangles;
ECStruct.opts.constants=ECStruct.Data.Cguess;
ECStruct.opts.ifit=1:ECStruct.Data.nsamp;
ECStruct.opts.iconst=1:length(ECStruct.Data.Cguess);

Cout=ECStruct.Data.Cguess;
ea=ECStruct.Data.eulerangles;
