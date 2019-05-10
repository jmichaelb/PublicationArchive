function [ECStruct,Cout,ea]=mkStrKspar(uncrt_prcnt)

ECStruct.Data.name='Kspar TestFitting';
ECStruct.Data.comments=['synthetic data'];
ECStruct.Data.nsamp=9; 
ECStruct.Data.sym='m';  % monoclinic example
ECStruct.Data.rho=2.600;  

ECStruct.Data.Cguess=[85
   50
   40
   -1
  160
   20
  -10
  165
   10
   20
  -12
   20
   30];

ECStruct.Data.Trust.constants=[
    50   100
    10    80
    10    80
   -20    20
   100   200
   -20    50
   -30    30
   100   200
   -30    30
    10    80
   -20    20
    10    80
    10    80
];

ECStruct.Data.Trust.eulerangles=[   10    10    10    10    10    10    10    10    10
                                    10    10    10    10    10    10    10    10    10
                                    40    40    40    40    40    40    40    40    40];


ECStruct.Data.eulerangles=[
  270   0     0     160    200   30    0  230    0
   90   90   55     125    125   90  155  150   25
  170  220   76      55    100  250  350  220   25];

% constrains from xrd data
ECStruct.Data.compliances=[ 8.11 3.27 3.69 0 0.2 0 0.01 0.01 0.01 .01 .01 .01]; 

% Al coating informations  - assumed isotropic 
atru = [1 0 0;0 1 0; 0 0 1]; 
rhou=2.70;
cmu=[  110.0   55.5   55.5     0     0     0
        55.5  110.0   55.5     0     0     0
        55.5   55.5  110.0     0     0     0
         0       0      0   27.0     0     0
         0       0      0      0   27.0    0
         0       0      0      0     0   27.0];
     
facSW=uncrt_prcnt/100;


%------------individual samples-------------------------------------------
%  all velocities below were calculated based on the nominal moduli
%  listed above.  Random error is then added to each velocity based 
%  on the uncert_prcnt specified in the input to this function
%    
% The mirror frequency (a velocity standard) specified below is used to calibrate
%   Impulsive Stimulated Light Scattering velocities and to determine the acoustic
%   wavelength which is needed to calculate the impact of the aluminum layer
%   on measured velocities.

ECStruct.Data.sample(1).name='Kspar 1 ';
ECStruct.Data.sample(1).mflg=0;
ECStruct.Data.sample(1).film.thickness=0.04;
ECStruct.Data.sample(1).film.density=rhou;
ECStruct.Data.sample(1).film.constants=cmu;
ECStruct.Data.sample(1).film.orientation=atru;
v_glass_surface=3.073; % the standard calibration velocity for this experiment

data1=[
         0    2.4369
         0    2.4369
   20.0000    2.6605
   30.0000    2.8140
   30.0000    2.8140
   50.0000    3.1614
  110.0000    2.8155
  130.0000    2.4797
  150.0000    2.1585
  170.0000    2.2567
  180.0000    2.4369
  160.0000    2.1465
  140.0000    2.2850
  120.0000    2.6829
   60.0000    3.2896
   40.0000    2.9928
   10.0000    2.5600
  105.0000    2.8012
  105.0000    2.8012
  115.0000    2.7657
  145.0000    2.2098
   90.0000    3.3837
   95.0000    3.3470
]; 

ECStruct.Data.sample(1).wavelength=v_glass_surface/1.830; %put in average mirror frequency
ECStruct.Data.sample(1).SWangles=data1(:,1);
ECStruct.Data.sample(1).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(1).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(1).SWivels=1:length(data1(:,2));


ECStruct.Data.sample(2).name='Kspar 2 ';
ECStruct.Data.sample(2).mflg=0;
ECStruct.Data.sample(2).film.thickness=0.04; 
ECStruct.Data.sample(2).film.density=rhou;
ECStruct.Data.sample(2).film.constants=cmu;
ECStruct.Data.sample(2).film.orientation=atru;
v_glass_surface=3.073; 
data1=[
          0    2.7102
   20.0000    2.8860
   20.0000    2.8860
   80.0000    2.8860
  100.0000    2.7100
  120.0000    2.5239
  160.0000    2.5240
  170.0000    2.5958
  150.0000    2.5229
  130.0000    2.5230
  110.0000    2.5957
   90.0000    2.8266
   30.0000    2.7472
   40.0000    2.3610
   20.0000    2.8860
   10.0000    2.8267
   35.0000    2.5640
  145.0000    2.5560
   60.0000    3.5469
   40.0000    3.5468
   35.0000    3.6136
   40.0000    3.5468
   30.0000    3.7549
   40.0000    3.5468
   55.0000    3.5240
   50.0000    3.5192
   60.0000    3.5469
   ];
 
ECStruct.Data.sample(2).wavelength=v_glass_surface/1.836; %put in average mirror frequency
ECStruct.Data.sample(2).SWangles=data1(:,1);
ECStruct.Data.sample(2).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));
ECStruct.Data.sample(2).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(2).SWivels=1:length(data1(:,2));%second way to exclude data via id



ECStruct.Data.sample(3).name='Kspar 3';
ECStruct.Data.sample(3).mflg=0;
ECStruct.Data.sample(3).film.thickness=0.04; %if needed change thickness here
ECStruct.Data.sample(3).film.density=rhou;
ECStruct.Data.sample(3).film.constants=cmu;
ECStruct.Data.sample(3).film.orientation=atru;
v_glass_surface=3.073; %measured 6/2/05 against water, believed accurate to 0.1%

data1=[ 
    20.0000    2.1598
   40.0000    2.4780
   60.0000    2.7166
   80.0000    2.7056
  140.0000    2.6970
  160.0000    2.6256
  170.0000    2.4374
  150.0000    2.7151
   90.0000    2.8260
   70.0000    2.6918
   50.0000    2.6548
   30.0000    2.2859
   95.0000    2.9271
  125.0000    2.7292
  130.0000    2.6955
  165.0000    2.5370
   65.0000    2.7064
   10.0000    2.1474
        ];

ECStruct.Data.sample(3).wavelength=v_glass_surface/1.836; %put in average mirror frequency
ECStruct.Data.sample(3).SWangles=data1(:,1);
ECStruct.Data.sample(3).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(3).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(3).SWivels=1:length(data1(:,2));%second way to exclude data via id


ECStruct.Data.sample(4).name='Kspar 4 ';
ECStruct.Data.sample(4).mflg=0;
ECStruct.Data.sample(4).film.thickness=0.04; %if needed change thickness here
ECStruct.Data.sample(4).film.density=rhou;
ECStruct.Data.sample(4).film.constants=cmu;
ECStruct.Data.sample(4).film.orientation=atru;
v_glass_surface=3.073;
data1=[
         0    2.5650
   20.0000    2.2407
   40.0000    2.2730
   60.0000    2.5950
   80.0000    2.7605
  100.0000    2.8286
  120.0000    3.0625
  140.0000    2.6791
  160.0000    2.7311
  180.0000    2.5650
  170.0000    2.7183
  180.0000    2.5650
  150.0000    2.6622
  130.0000    2.8752
  110.0000    2.9426
   90.0000    2.7715
   70.0000    2.7275
   50.0000    2.4221
   30.0000    2.2025
   10.0000    2.3775
  130.0000    3.3557];  

    
ECStruct.Data.sample(4).wavelength=v_glass_surface/1.837; %put in average mirror frequency
ECStruct.Data.sample(4).SWangles=data1(:,1);
ECStruct.Data.sample(4).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(4).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(4).SWivels=1:length(data1(:,2));%second way to exclude data via id


ECStruct.Data.sample(5).name='Kspar 5';
ECStruct.Data.sample(5).mflg=0;
ECStruct.Data.sample(5).film.thickness=0.04; %if needed change thickness here
ECStruct.Data.sample(5).film.density=rhou;
ECStruct.Data.sample(5).film.constants=cmu;
ECStruct.Data.sample(5).film.orientation=atru;
v_glass_surface=3.073; 

data1=[ 
         0    2.2077
   40.0000    2.7448
   60.0000    2.6518
   80.0000    3.0168
  100.0000    2.8795
  120.0000    2.7627
  120.0000    2.7627
  140.0000    2.6713
  160.0000    2.3406
  180.0000    2.2077
  170.0000    2.2251
  150.0000    2.5089
  130.0000    2.7539
  110.0000    2.7925
   90.0000    3.0109
   70.0000    2.7539
   50.0000    2.6955
   30.0000    2.6650
   10.0000    2.3200
   60.0000    3.4495
   70.0000    3.4072];

ECStruct.Data.sample(5).wavelength=v_glass_surface/1.837; %put in average mirror frequency
ECStruct.Data.sample(5).SWangles=data1(:,1);
ECStruct.Data.sample(5).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(5).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(5).SWivels=1:length(data1(:,2));%second way to exclude data via id


ECStruct.Data.sample(6).name='Kspar 6 ';
ECStruct.Data.sample(6).mflg=0;
ECStruct.Data.sample(6).film.thickness=0.04; %if needed change thickness here
ECStruct.Data.sample(6).film.density=rhou;
ECStruct.Data.sample(6).film.constants=cmu;
ECStruct.Data.sample(6).film.orientation=atru;
v_glass_surface=3.073; 
data1=[
         0    3.0434
   60.0000    2.7272
   80.0000    2.5598
  100.0000    2.5276
  140.0000    2.6411
  160.0000    2.9594
  180.0000    3.0434
  170.0000    3.0432
  150.0000    2.8302
  140.0000    2.6411
  135.0000    2.5173
  130.0000    2.4118
  125.0000    2.3606
  110.0000    2.5557
   90.0000    2.5199
   70.0000    2.6369
   10.0000    3.1139
   50.0000    2.7845
   50.0000    2.7845];      

ECStruct.Data.sample(6).wavelength=v_glass_surface/1.833; %put in average mirror frequency
ECStruct.Data.sample(6).SWangles=data1(:,1);
ECStruct.Data.sample(6).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(6).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(6).SWivels=1:length(data1(:,2));%second way to exclude data via id



ECStruct.Data.sample(7).name='Kspar 7';
ECStruct.Data.sample(7).mflg=0;
ECStruct.Data.sample(7).film.thickness=0.04; %if needed change thickness here
ECStruct.Data.sample(7).film.density=rhou;
ECStruct.Data.sample(7).film.constants=cmu;
ECStruct.Data.sample(7).film.orientation=atru;
v_glass_surface=3.073; 
data1=[
         0    2.9210
   20.0000    2.9235
   40.0000    2.7219
   60.0000    2.7224
   80.0000    2.4080
  140.0000    2.7225
  160.0000    2.7220
  180.0000    2.9190
  170.0000    2.7285
  150.0000    2.7373
  155.0000    2.7283
  130.0000    2.6020
   70.0000    2.6023
   50.0000    2.7373
  165.0000    2.7232
   10.0000    2.9579
   80.0000    3.4603
   95.0000    3.5948
   90.0000    2.2400
   90.0000    3.5632
  140.0000    3.4566
  100.0000    3.6040
  130.0000    3.4339];
   
ECStruct.Data.sample(7).wavelength=v_glass_surface/1.831; %put in average mirror frequency
ECStruct.Data.sample(7).SWangles=data1(:,1);
ECStruct.Data.sample(7).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(7).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(7).SWivels=1:length(data1(:,2));%second way to exclude data via id


ECStruct.Data.sample(8).name='Kspar 8 ';
ECStruct.Data.sample(8).mflg=0;
ECStruct.Data.sample(8).film.thickness=0.04; %if needed change thickness here
ECStruct.Data.sample(8).film.density=rhou;
ECStruct.Data.sample(8).film.constants=cmu;
ECStruct.Data.sample(8).film.orientation=atru;
v_glass_surface=3.073;
data1=[
         0    2.6827
   20.0000    2.6373
   40.0000    2.5242
   60.0000    2.5121
  100.0000    2.3310
  120.0000    2.4950
  140.0000    2.6717
  170.0000    2.6771
  130.0000    2.5800
  110.0000    2.4020
   90.0000    2.3063
   70.0000    2.4352
   50.0000    2.5376
   30.0000    2.5578
   15.0000    2.6690
    5.0000    2.6881
   25.0000    2.5965
   35.0000    2.5321
   75.0000    2.3928
   85.0000    2.3250
  180.0000    2.6827
  160.0000    2.7230
  150.0000    2.7420
];  

ECStruct.Data.sample(8).wavelength=v_glass_surface/1.834; %put in average mirror frequency
ECStruct.Data.sample(8).SWangles=data1(:,1);
ECStruct.Data.sample(8).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(8).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(8).SWivels=1:length(data1(:,2));%second way to exclude data via id


ECStruct.Data.sample(9).name='Kspar 9 ';
ECStruct.Data.sample(9).mflg=0;
ECStruct.Data.sample(9).film.thickness=0.04; %if needed change thickness here
ECStruct.Data.sample(9).film.density=rhou;
ECStruct.Data.sample(9).film.constants=cmu;
ECStruct.Data.sample(9).film.orientation=atru;
v_glass_surface=3.073; 

data1=[
         0    2.6271
   20.0000    2.5493
   40.0000    2.3264
   60.0000    2.1467
  100.0000    2.4530
  120.0000    2.5935
  140.0000    2.6851
  160.0000    2.7371
  180.0000    2.6271
  170.0000    2.6851
  170.0000    2.6851
  150.0000    2.7372
  130.0000    2.6271
  110.0000    2.5493
   90.0000    2.3266
   70.0000    2.1468
   50.0000    2.2130
   30.0000    2.4529
   10.0000    2.5935
  145.0000    2.7153 
];

ECStruct.Data.sample(9).wavelength=v_glass_surface/1.835; %put in average mirror frequency
ECStruct.Data.sample(9).SWangles=data1(:,1);
ECStruct.Data.sample(9).SWvelocities=data1(:,2).*(1+facSW*randn(size(data1(:,1))));;
ECStruct.Data.sample(9).SWsigvels=1.2*facSW*data1(:,2);
ECStruct.Data.sample(9).SWivels=1:length(data1(:,2));%second way to exclude data via id
%_____________________ End of surface wave data entry

% __________ define variables in the opts section
ECStruct.opts.ifit=1:ECStruct.Data.nsamp; % default to fit all data
ECStruct.opts.iconst=1:length(ECStruct.Data.Cguess); % default to optimize all moduli
ECStruct.opts.vrange=.2; % search region about each surface wave measurement
ECStruct.opts.mflg=0; % 0 not a mirror plane, 1 if mirror plane
ECStruct.opts.vmin=2; % define upper and lower plausible surface wave velocities
ECStruct.opts.vmax=4;
ECStruct.opts.constrflg='y'; % y if you have x ray compliances
ECStruct.opts.eulerangles=ECStruct.Data.eulerangles; % put current "best" Euler angles in the opts section


Cout=ECStruct.Data.Cguess;
ea=ECStruct.Data.eulerangles;



