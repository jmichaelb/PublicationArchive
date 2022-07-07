function [ECStruct, Cout, ea]=mkStrFoMao2p7_500K
% data from:
%Elasticity of single-crystal olivine at high pressures and temperatures
%Zhu Mao, DaweiFanc, Jung-Fu Lin, JingYang Sergey N.Tkachev, KirillZhuravlev, Vitali B.Prakapenka
% Earth and Planetary Science Letters 426(2015)204–215
% excel spreadsheet of velocity data provided by author.  No euler angles
% were provided in spreadsheet or publication.  The euler angles given below were found by fitting.

ECStruct.Data.name='Mao et al 2.7 GPa 500 K';
ECStruct.Data.dcosflg=0;  % samples are organized as data in two planes.
ECStruct.Data.nsamp=2;
ECStruct.Data.sym='o';
ECStruct.Data.rho=3.381;

% the Cij reported in the publication
ECStruct.Data.Cguess=[ 332   77  79   204	 82    240	 66.1	 79.4	  82 ];% [ 340   80   80 200  80   238    70   80     80];
%                   [ 330   76   79  204.7  84   241    65.9    79.     80];
                     

% estimated bounds on moduli
ECStruct.Data.Trust.constants=[  300   400
                       60   100
                      60  100
                     180  300
                       60   100
                       200 300
                      50 100
                       50 100
                      50 100
];
ECStruct.Data.Trust.eulerangles=[180 ;180; 180]*ones(1,2);



% the following euler angles were found to fit the Cij reported 
ECStruct.Data.eulerangles=[   172.1002  -61.1513
                              64.8894  -97.6554
                              91.6518  146.5819];

% here is other set that gives the same misfit and are essentially symmetry related.
%    -7.7908  118.7388
%   -64.8934   96.5934
%   -88.4385  146.0553

%Sample 1
ECStruct.Data.sample(1).name='Platelet 1';  
% data were taken at rotations of 10 degrees:
ECStruct.Data.sample(1).BWangles= [0:10:170]';

% the following data were copied from the spreadsheet.

%		azimuthal angle	velocity	error	density
%plane 1	Vp	
Vp_p1= [ 0	7.729	0.025	3.381
		10	7.74	0.025	3.381
		20	7.795	0.031	3.381
		30	7.918	0.034	3.381
		40	7.965	0.029	3.381
		50	8.131	0.033	3.381
		60	8.253	0.026	3.381
		70	8.368	0.034	3.381
		80	8.422	0.029	3.381
		90	8.438	0.051	3.381
		100	8.483	0.032	3.381
		110	8.424	0.029	3.381
		120	8.309	0.036	3.381
		130	8.176	0.028	3.381
        140 NaN      NaN     NaN 
        150 NaN      NaN     NaN 
		160	7.962	0.031	3.381
		170	7.848	0.026	3.381];

	Vs1_p1=[
    	0	4.423	0.015	3.381
		10	4.427	0.015	3.381
		20	4.482	0.016	3.381
		30	4.514	0.017	3.381
		40	4.535	0.015	3.381
		50	4.544	0.017	3.381
		60	4.533	0.016	3.381
		70	4.516	0.014	3.381
		80	4.507	0.014	3.381
		90	4.49	0.014	3.381
		100	4.541	0.015	3.381
		110	4.577	0.015	3.381
		120	4.586	0.015	3.381
		130	4.569	0.016	3.381
		140	4.598	0.017	3.381
		150	4.561	0.015	3.381
		160	4.514	0.015	3.381
		170	4.468	0.014	3.381];

	Vs2_p1=[
    	0	4.876	0.035	3.381
		10	4.905	0.038	3.381
        20  NaN      NaN     NaN 
        30 NaN      NaN     NaN 
        40 NaN      NaN     NaN 
        50 NaN      NaN     NaN 
        60 NaN      NaN     NaN
        70 NaN      NaN     NaN
        80 NaN      NaN     NaN 
        90 NaN      NaN     NaN 
        100 NaN      NaN     NaN 
        110 NaN  NaN  NaN
        120 NaN      NaN     NaN 
        130 NaN      NaN     NaN 
        140 NaN      NaN     NaN 
        150 NaN      NaN     NaN 
        160 NaN      NaN     NaN 
		170	4.961	0.021	3.381];


    % the following two lines organizes the data in the format used by
    % Velocities2Cij
    % Bodywave velocities listed in ascending order - three polarizations for
    % each angle - "NaN" (not a number) denotes an unobserved mode.
    ECStruct.Data.sample(1).BWvelocities=[Vs1_p1(:,2) Vs2_p1(:,2) Vp_p1(:,2) ];
    ECStruct.Data.sample(1).BWuncertainties=[Vs1_p1(:,3) Vs2_p1(:,3) Vp_p1(:,3) ];


    %plane 2	
Vp_p2=[
         0	8.538	0.041	3.381
		10	8.54	0.043	3.381
		20	8.438	0.029	3.381
		30	8.441	0.031	3.381
        40 NaN  NaN   NaN
		50	8.39	0.037	3.381
		60	8.345	0.032	3.381
		70	8.42	0.030	3.381
		80	8.569	0.034	3.381
		90	8.797	0.034	3.381
		100	9.068	0.049	3.381
		110	9.229	0.054	3.381
		120	9.334	0.036	3.381
		130	9.368	0.037	3.381
		140	9.314	0.032	3.381
		150	9.169	0.030	3.381
		160	8.959	0.036	3.381
		170	8.718	0.042	3.381];

	Vs2_p2=[
    	0	5.405	0.018	3.381
		10	5.247	0.017	3.381
		20	5.034	0.017	3.381
		30	4.877	0.016	3.381
		40	4.858	0.019	3.381
		50	4.965	0.017	3.381
		60	5.124	0.017	3.381
		70	5.284	0.017	3.381
		80	5.347	0.016	3.381
		90	5.312	0.017	3.381
		100	5.198	0.020	3.381
		110	5.103	0.018	3.381
        120 NaN  NaN   NaN
		130	5.058	0.017	3.381
		140	5.103	0.022	3.381
		150	5.214	0.024	3.381
		160	5.354	0.017	3.381
		170	5.429	0.018	3.381];


	Vs1_p2=[
    	0	4.582	0.015	3.381
		10	4.47	0.020	3.381
		20	4.447	0.021	3.381
		30	4.454	0.015	3.381
		40	4.424	0.019	3.381
		50	4.504	0.034	3.381
		60	4.478	0.021	3.381
		70	4.521	0.019	3.381
        80 NaN  NaN   NaN
		90	4.69	0.024	3.381
		100	4.732	0.019	3.381
		110	4.766	0.018	3.381
		120	4.728	0.016	3.381
		130	4.748	0.016	3.381
		140	4.758	0.015	3.381
		150	4.737	0.015	3.381
		160	4.715	0.016	3.381
		170	4.635	0.016	3.381];
ECStruct.Data.sample(2).name='Platelet 2';  
% Bodywave velocities listed in ascending order - three polarizations for
% each angle - "NaN" (not a number) denotes an unobserved mode.
% the following two lines organizes the data in the format used by
% Velocities2Cij
    ECStruct.Data.sample(2).BWvelocities=[Vs1_p2(:,2) Vs2_p2(:,2) Vp_p2(:,2) ];
    ECStruct.Data.sample(2).BWuncertainties=[Vs1_p2(:,3) Vs2_p2(:,3) Vp_p2(:,3) ];

% data were taken at rotations of 10 degrees:
ECStruct.Data.sample(2).BWangles= [0:10:170]';


  
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
