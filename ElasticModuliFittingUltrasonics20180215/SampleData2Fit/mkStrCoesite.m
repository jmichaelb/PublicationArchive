function [ECStruct,Co,ea]=mkStrCoesite(Cflg)
% Weidner and Carlton (1977) data copied from paper

ECStruct.Data.name='Coesite';

% all data reported with direction cosines
ECStruct.Data.dcosflg=1;

% group all data together
ECStruct.Data.nsamp=1;

% set symmetry and density
ECStruct.Data.sym='m';
ECStruct.Data.rho=2.92;

% starting estimates for the moduli (cyclic order)
if strcmp(Cflg,'p')
ECStruct.Data.Cguess=[160.8 82.1 102.9 -36.2 230.4 35.6 2.6 231.6 -39.3 67.8 9.9 73.3 58.8]';
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
ECStruct.Data.Trust.constants=[  50   300
                       10   100
                       10   150
                      -60   100
                       50   400
                       10   150
                      -60    60
                       50   400
                      -50    50
                       15   100
                      -50    50
                       10   100
                       10   100
];

% data copied directly from paper
%     direction cosines  Vp Vs1 Vs2
data=[-.208 .978 .003 8.76 6.41 4.49
.999 .035 -.014 7.69 4.58 nan
-.114 .991 .066 8.73 nan nan
.114 .992 -.062 9.00 nan nan
.334 .924  -.185 8.62 nan 4.49
.531 .794 -.296 8.22 nan nan
.692 .610 -.387 8.51 6.35 3.98
.806 .384 -.451 8.95 nan nan
.864 .131 -.484 9.51 nan nan
.279 .959 .069 8.62 6.19 nan
.502 .851 .153 8.62 4.53 nan
.691 .686 .226 8.29 nan nan
.833 .475 .284 7.51 nan nan
.918 .231 .332 7.01 4.53 nan
.940 -.029 .339 6.77 4.78 nan
-.156 -.988 0 8.85 nan nan
.199 -.980 0 8.69 nan nan
.367 -.930 0 8.61 5.96 4.19
.523 -.853 0 8.37 nan nan
.595 -.804 0 8.34 nan 4.56
.663 -.749 0 8.18 5.23 3.86
.783 -.623 0 8.05 5.39 3.94
.875 -.485 0 7.89 5.28 3.70
.948 -.317 0 7.73 4.99 4.02
.989 -.148 0 7.70 nan 4.35
1 .026 0 7.65 nan 4.35
.980 .199 0 7.73 5.15 3.86
.930 .367 0 7.73 5.35 3.83
.853 .522 0 7.97 5.51 3.83
.804 .595 0 8.21 4.91 3.54
.749 .663 0 8.34 5.07 4.02
.623 .783 0 8.45 nan nan
.477 .879 0 8.61 6.04 3.99
.326 .946 0 8.61 6.60 5.15
.156 .988 0 8.94 nan 3.94
-.964 0 .266 9.16 nan nan
-.903 0 .429 9.39 nan nan
-.815 0 .579 9.89 nan nan
-.702 0 .712 10.19 nan nan
-.335 0 .942 9.76 nan nan
-.167 0 .986 9.60 6.00 nan
.007 0 1 8.89 nan nan
.181 0 .984 8.54 nan 4.54
.429 0 .903 7.99 nan 4.35
.648 0 .762 7.32 5.34 nan
.823 0 .568 6.80 4.81 nan
.942 0 .335 6.91 4.63 nan
.997 0 .080 7.53 4.60 nan
.984 0 -.181 8.45 nan nan
-.730 -.227 .644 9.76 7.35 nan
-.294 -.380 .877 9.25 nan nan
-.037 -.420 .907 8.66 nan nan
.221 -.431 .875 8.20 nan 5.25
.465 -.413 .783 7.76 nan nan
.677 -.367 .638 7.49 nan nan
.843 -.295 .450 7.21 nan nan
.952 -.204 .230 7.24 nan 4.40
.971 .014 -.239 8.73 nan nan
.120 .797 .592 8.53 nan 4.67
.445 .878 .175 8.73 4.53 nan
-.217 .514 .830 8.91 nan 4.64
-.191 .973 .129 8.95 nan 4.64
-.318 .885 .339 8.67 5.86 nan
-.477 .610 .633 8.78 6.08 nan
-.549 .134 .825 10.39 nan nan 
-.473 -.377 .796 9.88 nan nan
-.385 -.603 .699 9.02 nan nan
-.271 -.787 .554 8.31 nan 4.64
-.138 -.918 .372 8.67 nan 4.80
.004 -.987 .163 8.95 4.86 nan
.278 -.922 -.271 8.47 nan 4.64
.136 -.971 .195 8.91 6.19 4.38
-.057 -.773 .632 8.36 nan 4.64
-.234 -.367 .900 9.51 nan 4.49
-.559 .626 .544 9.22 6.63 5.86
-.708 .152 .690 9.91 nan nan
-.668 -.362 .650 9.79 nan nan
-.578 -.591 .563 9.09 nan nan
-.448 -.780 .437 8.36 nan nan
-.109 -.988 .106 9.07 nan nan
.260 -.932 -.253 8.62 nan 4.91
];


  ECStruct.Data.sample(1).name='Coesite';
%set direction cosines:
    ECStruct.Data.sample(1).dcos= data(:,1:3);
    
% here velocities are listed in ascending order
  ECStruct.Data.sample(1).BWvelocities= data(:,[6 5 4]);
  ECStruct.Data.sample(1).BWangles= []; % angles are not used here
  
% the overall misfit is 0.15 km/s.  However, a check of shear waves and
% compressional waves indicate a slightly larger scatter for the shear
% waves. Thus, uncertainties are assigned as 0.18 for the shear and 0.13
% for the compressional waves
  ECStruct.Data.sample(1).BWuncertainties=[.18*ones(81,2) .13*ones(81,1)];

  
% set necessary options:
ECStruct.opts.constants=ECStruct.Data.Cguess;
ECStruct.opts.ifit=1:ECStruct.Data.nsamp;
ECStruct.opts.iconst=1:length(ECStruct.Data.Cguess);

Co=ECStruct.Data.Cguess;
ea=[];

% - according to the paper all data deviating from the model by more
% than 0.5 km/s were excluded from the fit.  So to determine which points
% are to be excluded, the data as defined above were loaded and using
% published moduli, deviations were calculated and pasted in below along
% with the related direction cosines. Below this matrix, a set of steps
% identify the data to exclude and reform the output data set.

%         direction cosines                deviations
r=[ -0.2080    0.9780    0.0030    0.1751    1.4480   -0.0757
    0.9990    0.0350   -0.0140       NaN   -0.0032   -0.0818
    %
   -0.1140    0.9910    0.0660       NaN       NaN   -0.1195
    0.1140    0.9920   -0.0620       NaN       NaN    0.1425
    0.3340    0.9240   -0.1850   -0.1238       NaN   -0.0240
    0.5310    0.7940   -0.2960       NaN       NaN   -0.2222
    0.6920    0.6100   -0.3870    0.0561    0.5319   -0.1750
    0.8060    0.3840   -0.4510       NaN       NaN   -0.2689
    0.8640    0.1310   -0.4840       NaN       NaN   -0.0569
    %
    0.2790    0.9590    0.0690       NaN    1.2897   -0.2033
    0.5020    0.8510    0.1530       NaN   -0.1550    0.0113
    0.6910    0.6860    0.2260       NaN       NaN    0.0683
    0.8330    0.4750    0.2840       NaN       NaN   -0.1553
    0.9180    0.2310    0.3320       NaN   -0.0686   -0.0532
    0.9400   -0.0290    0.3390       NaN    0.0259    0.0078
    %
   -0.1560   -0.9880         0       NaN       NaN   -0.0097
    0.1990   -0.9800         0       NaN       NaN   -0.1514
    0.3670   -0.9300         0   -0.2160    1.1075   -0.1246
    0.5230   -0.8530         0       NaN       NaN   -0.2054
    0.5950   -0.8040         0    0.3691       NaN   -0.1295
    0.6630   -0.7490         0   -0.2048    0.1000   -0.1771
    0.7830   -0.6230         0    0.0784    0.0973   -0.0756
    0.8750   -0.4850         0   -0.0940   -0.0402   -0.0383
    0.9480   -0.3170         0    0.1097   -0.1660   -0.0520
    0.9890   -0.1480         0    0.1554       NaN   -0.0278
    1.0000    0.0260         0   -0.1021       NaN   -0.0687
    0.9800    0.1990         0   -0.2369    0.1998   -0.0088
    0.9300    0.3670         0   -0.0253    0.1259   -0.0847
    0.8530    0.5220         0    0.0353    0.1842   -0.0009
    0.8040    0.5950         0   -0.2935   -0.3989    0.1326
    0.7490    0.6630         0    0.1091   -0.1833    0.1472
    0.6230    0.7830         0       NaN       NaN    0.0215
    0.4770    0.8790         0   -0.3767    1.1632   -0.0176
    0.3260    0.9460         0    0.7580    1.7215   -0.1640
    0.1560    0.9880         0   -0.3413       NaN    0.0803
    %
   -0.9640         0    0.2660       NaN       NaN    0.3347
   -0.9030         0    0.4290       NaN       NaN   -0.0394
   -0.8150         0    0.5790       NaN       NaN    0.0131
   -0.7020         0    0.7120       NaN       NaN    0.0517
   -0.3350         0    0.9420       NaN       NaN   -0.1663
   -0.1670         0    0.9860       NaN    1.3076    0.0516
    0.0070         0    1.0000       NaN       NaN   -0.1710
    0.1810         0    0.9840   -0.1025       NaN   -0.0021
    0.4290         0    0.9030    0.0530       NaN    0.1541
    0.6480         0    0.7620       NaN    0.3115    0.0138
    0.8230         0    0.5680       NaN   -0.1195   -0.0915
    0.9420         0    0.3350       NaN   -0.1255    0.1490
    0.9970         0    0.0800       NaN   -0.0344    0.1302
    0.9840         0   -0.1810       NaN       NaN   -0.0311
    %
   -0.7300   -0.2270    0.6440       NaN    2.7250   -0.1266
   -0.2940   -0.3800    0.8770       NaN       NaN   -0.1847
   -0.0370   -0.4200    0.9070       NaN       NaN   -0.1275
    0.2210   -0.4310    0.8750    0.9966       NaN   -0.0001
    0.4650   -0.4130    0.7830       NaN       NaN   -0.0787
    0.6770   -0.3670    0.6380       NaN       NaN   -0.0572
    0.8430   -0.2950    0.4500       NaN       NaN   -0.0007
    0.9520   -0.2040    0.2300    0.0011       NaN    0.1444
    0.9710    0.0140   -0.2390       NaN       NaN    0.0134
    %
    0.1200    0.7970    0.5920    0.2669       NaN    0.1510
    0.4450    0.8780    0.1750       NaN   -0.1951    0.0585
    %
   -0.2170    0.5140    0.8300    0.1089       NaN   -0.0585
   %
   -0.1910    0.9730    0.1290    0.1415       NaN    0.1649
   -0.3180    0.8850    0.3390       NaN    0.2585    0.1638
   -0.4770    0.6100    0.6330       NaN    0.2552   -0.1734
   -0.5490    0.1340    0.8250       NaN       NaN    0.2574
   -0.4730   -0.3770    0.7960       NaN       NaN    0.2018
   -0.3850   -0.6030    0.6990       NaN       NaN    0.0845
   -0.2710   -0.7870    0.5540   -0.0124       NaN   -0.0974
   -0.1380   -0.9180    0.3720    0.0973       NaN    0.0943
    0.0040   -0.9870    0.1630       NaN   -0.2027    0.1112

    0.2780   -0.9220   -0.2710   -0.0420       NaN   -0.1378
    0.1360   -0.9710    0.1950    0.1667    1.1731    0.0857
   -0.0570   -0.7730    0.6320    0.0235       NaN    0.0263
   -0.2340   -0.3670    0.9000    0.0091       NaN    0.1691

   -0.5590    0.6260    0.5440    1.6511    0.7638    0.3564
   -0.7080    0.1520    0.6900       NaN       NaN   -0.1352
   -0.6680   -0.3620    0.6500       NaN       NaN    0.1101
   -0.5780   -0.5910    0.5630       NaN       NaN    0.1048
   -0.4480   -0.7800    0.4370       NaN       NaN   -0.1016
   -0.1090   -0.9880    0.1060       NaN       NaN    0.2376
    0.2600   -0.9320   -0.2530    0.2302       NaN   -0.0159];


v=data(:,[6 5 4]);
dv=r(:,4:6);
idfnt=find(isfinite(v));
idb=find(abs(dv(idfnt))>.5);
v(idfnt(idb))=nan*ones(size(idb));
dv(idfnt(idb))=nan*ones(size(idb));
v=reshape(v,81,3);
dv=reshape(dv,81,3);


ECStruct.Data.sample(1).BWvelocities= v(:,[1 2 3]);


 