function [ECStruct, Cout, ea]=mkStrCPX(Cflg)
% Collins and Brown Elasticity of an upper mantle clinopyroxene, Phys Chem
% Minerals 26, 7-13, 1998

ECStruct.Data.name='Collins and Brown CPX 1998';
ECStruct.Data.dcosflg=0;  % samples are organized as data in three planes.
ECStruct.Data.nsamp=3;
ECStruct.Data.sym='m';
ECStruct.Data.rho=3.327;
if strcmp(Cflg,'p')
ECStruct.Data.Cguess=[ 237.8  83.5  80   9. 183.6  59.9   9.5 229.5  48.0  76.6   8.4  73  81.6 ];
else
ECStruct.Data.Cguess=[  100 50 50 0 100 50 0 100 0 50 0 50 50];
end
ECStruct.Data.Trust.constants=[  50   300
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
ECStruct.Data.Trust.eulerangles=ones(3,1)*[180 180 180];

ECStruct.Data.eulerangles=[    7.9  269.0  345.2
                              89.6   85.4    7.3
                               3.2  193.4  345.5];

fac=.004; % an uncertainty of 0.4% is associated with these measured velocities 
% The experimental results for three crystals are listed:

%Sample 1
ECStruct.Data.sample(1).name='a-normal cut Nov 30, 1996';  
% Bodywave velocities listed in ascending order - three polarizations for
% each angle - "NaN" (not a number) denotes an unobserved mode.
ECStruct.Data.sample(1).BWvelocities= [    
    4.1940       NaN    8.5540
    4.1790       NaN       NaN
    4.1690    4.8330    8.5330
    4.1940    4.8660    8.4330
       NaN    4.9290    8.2960
       NaN    4.9210    8.0890
       NaN    4.9840    7.9140
       NaN    5.0000    7.7020
       NaN    5.0510    7.5420
       NaN       NaN    7.4890
       NaN       NaN    7.4950
    4.3590    5.1460    7.5310
    4.2960    5.1570    7.6610
    4.2330    5.1520    7.7870
       NaN    5.0820    8.0330
    4.2590    4.9370    8.1650
    4.3140    4.8770    8.3610
    4.2700       NaN    8.4580
    4.2250       NaN    8.5450];

% data were taken at rotations of 10 degrees:
ECStruct.Data.sample(1).BWangles= [0:10:180]';

% uncertainties calculated and loaded into the structure:
ECStruct.Data.sample(1).BWuncertainties=fac*ECStruct.Data.sample(1).BWvelocities;
 

% Repeat steps for second sample:
ECStruct.Data.sample(2).name='b-normal cut July 11, 1997';
ECStruct.Data.sample(2).BWvelocities=[
       NaN    4.8540    7.7740
       NaN       NaN    7.1780
       NaN    5.2290    7.0060
       NaN    4.9780    7.1060
    4.4600       NaN    7.4410
    3.9300       NaN       NaN
    4.0510       NaN    8.2130
    4.3790       NaN    8.3670
    4.8060       NaN    8.5140
       NaN    5.1620       NaN
       NaN    5.2760    8.8150
       NaN    5.1570    9.0580
    4.8200       NaN       NaN
    4.4950       NaN    9.4140
    4.1300       NaN    9.4000
    3.9820       NaN       NaN
    4.0380       NaN       NaN
    4.3370       NaN    8.3300
       NaN    4.8190    7.7770];   
ECStruct.Data.sample(2).BWangles= [
    0
    10
    20
    25
    35
    50
    60
    70
    80
    90
   100
   110
   120
   130
   140
   150
   160
   170
   180];
ECStruct.Data.sample(2).BWuncertainties=fac*ECStruct.Data.sample(2).BWvelocities;
 
% Repeat steps for third sample:
ECStruct.Data.sample(3).name='c-normal cut June 17, 1997';
ECStruct.Data.sample(3).BWvelocities= [
    4.3960       NaN    8.3360
    4.3780       NaN    8.3450
    4.3770       NaN    8.3770
       NaN       NaN    8.3550
    4.3620       NaN    8.4230
    4.3860       NaN    8.4280
    4.3850       NaN    8.3720
    4.3970       NaN    8.3540
    4.3200    4.6930    8.1900
    4.4600       NaN    7.9690
    4.5600       NaN    7.7540
       NaN       NaN    7.5390
       NaN       NaN       NaN
       NaN       NaN    7.5300
    4.4670       NaN    7.7460
    4.4330       NaN    7.9410
    4.3260    4.5940    8.1590
    4.3120       NaN    8.2890];

  ECStruct.Data.sample(3).BWangles= [
           0
    10
    20
    30
    40
    50
    60
    70
    80
    90
   100
   110
   120
   130
   140
   150
   160
   170];
  ECStruct.Data.sample(3).BWuncertainties=fac*ECStruct.Data.sample(3).BWvelocities;

  
% A few things are loaded into the .opts side of the structure. These
% include: "ifit" which is a list of indexes of sampels to fit and "iconst"
% which is a list of moduli to fit.  Here we are fitting all data for all
% moduli. The moduli and euler angles loaded into the .opts side are
% modified during fitting as specified by the parameters used in the call
% to Velocities2Cij.  The values of moduli and euler angles loaded into the
% .Data side are not modified during fitting to retain a record of where
% the optimization started.
ECStruct.opts.constants=ECStruct.Data.Cguess;
ECStruct.opts.ifit=1:ECStruct.Data.nsamp;
ECStruct.opts.iconst=1:length(ECStruct.Data.Cguess);
ECStruct.opts.eulerangles=ECStruct.Data.eulerangles;

Cout=ECStruct.Data.Cguess;
ea=ECStruct.Data.eulerangles;
