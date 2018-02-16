function [ECStruct,Cout,ea]=mkStrHornblende(Cflg)
% data from Brown and Abramson "Elasticity of Amphiboles" submitted

ECStruct.Data.name='Example Hornblende Data ';
% Chemical Name: Potassian-Magnesio-Hastingitic Hornblende

% number of samples (slices with different orientations)
ECStruct.Data.nsamp=3;

%symmetry of crystal
ECStruct.Data.sym='m';

%density in gm/cc
ECStruct.Data.rho=3.2611;

% either set inital moduli to a previously determined set or to default
% values
if strcmp(Cflg,'p')
ECStruct.Data.Cguess=[ 132.8  
        53.6  
   48.1  
   -.9  
  189.2   
   61.1  
   -8.4  
  227.4  
     -30.9  
   73.4   
   4.3
   47.3  
   48.7  ];  
else
 ECStruct.Data.Cguess=[ 100  
                         50  
                         50  
                          0  
                        200   
                         50 
                          0  
                        200 
                          0  
                         50   
                          0
                         50 
                         50  ];  
end

% modify trust region for specific examples
ECStruct.Data.Trust.constants=[   100   150
                                   30    80
                                   30    80
                                  -20    40
                                  150   250
                                   30    80
                                  -20    50
                                  200   300
                                  -50    50
                                   30    80
                                  -30    30
                                  -40    80
                                  -20    80];
                              
% euler angles                              
ECStruct.Data.eulerangles=[
  180.54  274.03  359.62
   90.81   90.60   16.12
  260.34  230.46  156.73];

ECStruct.Data.Trust.eulerangles=[20;20;20]*ones(1,ECStruct.Data.nsamp);

% assume percentage uncertainty for velocities (otherwise modify line below
% where uncertainties are entered
fac=.0025;  %  percent uncertainty in velocities - experiment dependent

% bodywave velocities:
% nominal a normal sample
ECStruct.Data.sample(1).name='Hornblende   a* normal';
ECStruct.Data.sample(1).BWvelocities=[
       NaN       NaN    7.6620
       NaN       NaN    7.6190
       NaN       NaN    7.6740
       NaN    4.7410    7.7460
       NaN    4.7660    7.8780
       NaN       NaN    7.9770
       NaN    4.8070    8.1230
       NaN    4.8030       NaN
       NaN       NaN    8.3560
       NaN       NaN    8.4340
       NaN       NaN    8.4840
       NaN       NaN       NaN
       NaN    4.7780    8.3750
    3.6330    4.7870    8.2660
       NaN    4.7870    8.1230
       NaN    4.7880    7.9740
       NaN       NaN    7.7910
       NaN       NaN    7.6610
       NaN       NaN    7.6590];
ECStruct.Data.sample(1).BWangles=0:10:180;
%use constant percentage uncertainty (can be modified)
ECStruct.Data.sample(1).BWuncertainties=fac*ECStruct.Data.sample(1).BWvelocities;

% nominal b normal sample
ECStruct.Data.sample(2).name='Hornblende b normal';
ECStruct.Data.sample(2).BWvelocities= [
         NaN       NaN    7.6520
       NaN    4.5100    7.1120
       NaN    4.5130    6.6860
       NaN       NaN    6.4540
       NaN       NaN    6.4120
       NaN       NaN    6.3890
       NaN       NaN    6.3530
       NaN       NaN    6.2730
       NaN       NaN    6.2590
       NaN       NaN    6.4990
       NaN       NaN    6.9960
    4.1400       NaN    7.5770
    3.8290       NaN    8.0870
    3.5730       NaN    8.4560
    3.4660       NaN    8.6630
       NaN       NaN       NaN
    3.6840       NaN       NaN
       NaN       NaN       NaN
       NaN       NaN    7.6500];    
ECStruct.Data.sample(2).BWangles= 0:10:180;
%use constant percentage uncertainty (can be modified)
ECStruct.Data.sample(2).BWuncertainties=fac*ECStruct.Data.sample(2).BWvelocities;

% nominal c normal sample
ECStruct.Data.sample(3).name='Hornblende c normal';
ECStruct.Data.sample(3).BWvelocities=[ 
       NaN       NaN       NaN
       NaN    4.4360    6.6090
       NaN       NaN       NaN
       NaN    4.4600    6.5500
       NaN       NaN    6.6010
       NaN       NaN    6.7440
       NaN       NaN    6.8930
    3.9330    4.6340    7.0670
       NaN       NaN    7.2580
       NaN       NaN    7.4430
       NaN       NaN    7.5400
       NaN       NaN       NaN
       NaN       NaN       NaN
       NaN       NaN    7.5360
       NaN       NaN    7.3710
       NaN       NaN    7.2060
    3.9500    4.5950    6.9970
       NaN    4.5480    6.8320
       NaN       NaN    6.7030];    
 
ECStruct.Data.sample(3).BWangles= 0:10:180;
%use constant percentage uncertainty (can be modified)
ECStruct.Data.sample(3).BWuncertainties=fac*ECStruct.Data.sample(3).BWvelocities;
 
% Surfacewave specific below here
%aluminum coating properties
% thickness of coating in microns
thick=.04; 
% orientation matrix for coating (since the coating is isotropic, this is
% just an identity matrix)
atru = [1 0 0;0 1 0; 0 0 1]; 
%density of coating in gm/cc
rhou=2.70;
% isotropic moduli for aluminum (GPa)
cmu=[  110.0   55.5   55.5     0     0     0
        55.5  110.0   55.5     0     0     0
        55.5   55.5  110.0     0     0     0
         0       0      0   27.0     0     0
         0       0      0      0   27.0    0
         0       0      0      0     0   27.0];

% glass standard used to determine experimental wavelength   
v_glass_surface=3.073; %measured 6/2/05 against water, believed accurate to 0.1%

% Surface wave veolcities:
facSW=.0025;
% nominal a normal sample
ECStruct.Data.sample(1).wavelength=v_glass_surface/1.2112;
ECStruct.Data.sample(1).SWvelocities=[
    3.5670
    3.5640
    3.5710
    3.5840
    3.5700
    3.5350
    3.4860
    3.4360
    3.4020
    3.3920
    3.3710
    3.3860
    3.3940
    3.4400
    3.4850
    3.5400
    3.5720
    3.5750
    3.5700];

ECStruct.Data.sample(1).SWangles=(0:10:180)';
%use constant percentage uncertainty (can be modified)
ECStruct.Data.sample(1).SWsigvels=facSW*ECStruct.Data.sample(1).SWvelocities;

ECStruct.Data.sample(1).SWivels=1:length(ECStruct.Data.sample(1).SWvelocities);
ECStruct.Data.sample(1).mflg=0;
ECStruct.Data.sample(1).film.thickness=thick; 
ECStruct.Data.sample(1).film.density=rhou;
ECStruct.Data.sample(1).film.constants=cmu;
ECStruct.Data.sample(1).film.orientation=atru;

% nominal b normal sample
ECStruct.Data.sample(2).wavelength=v_glass_surface/1.211;
ECStruct.Data.sample(2).SWvelocities=[
    3.7950
    3.6550
    3.5990
    3.5580
    3.6120
    3.6350
    3.7390
    3.8470
    3.9590
    4.0360
    4.0680
    4.3180
    4.2300
    4.1140
    3.7590];    
ECStruct.Data.sample(2).SWangles=[0:10:30 50:10:110 140:10:160 180]';
ECStruct.Data.sample(2).SWsigvels=facSW*ECStruct.Data.sample(2).SWvelocities;

% ivels is list of velocities to fit - remove indexes to exclude some data
% from being fit
ECStruct.Data.sample(2).SWivels=1:length(ECStruct.Data.sample(2).SWvelocities);
ECStruct.Data.sample(2).mflg=0;
ECStruct.Data.sample(2).film.thickness=thick;  
ECStruct.Data.sample(2).film.density=rhou;
ECStruct.Data.sample(2).film.constants=cmu;
ECStruct.Data.sample(2).film.orientation=atru;

% nominal c normal sample
ECStruct.Data.sample(3).wavelength=v_glass_surface/1.208;
ECStruct.Data.sample(3).SWvelocities=[
    3.7290
    3.6470
    3.6360
    3.6530
    3.6700
    3.7420
    3.7920
    3.8570
    4.1520
    4.2330
    4.2570
    4.2440
    4.2080
    3.8600
    4.1100
    3.7800
    3.7170];    
ECStruct.Data.sample(3).SWangles=[0:10:60 90 90:10:130 150 150 170 180]';
ECStruct.Data.sample(3).SWsigvels=facSW*ECStruct.Data.sample(3).SWvelocities;
ECStruct.Data.sample(3).SWivels=1:length(ECStruct.Data.sample(3).SWvelocities);
ECStruct.Data.sample(3).mflg=0;
ECStruct.Data.sample(3).film.thickness=thick;
ECStruct.Data.sample(3).film.density=rhou;
ECStruct.Data.sample(3).film.constants=cmu;
ECStruct.Data.sample(3).film.orientation=atru;


% set options
% samples to fit
ECStruct.opts.ifit=[1 2 3];
% moduli to fit
ECStruct.opts.iconst=[1 2 3 4 5 6 7 8 9 10 11 12 13];

% trial values for euler angles and moduli
ECStruct.opts.eulerangles=ECStruct.Data.eulerangles;
ECStruct.opts.constants=ECStruct.Data.Cguess;

% range of angles for plots
ECStruct.opts.pltRange=180;

% if axes compresibilities are used as constraints set to 'y'
ECStruct.opts.constrflg='n';

% for surface waves set the search range:
ECStruct.opts.vrange=0.2;
ECStruct.opts.vmin=3.3;
ECStruct.opts.vmax=4.5;

% return starting values for moduli and euler angles
Cout=ECStruct.Data.Cguess;
ea=ECStruct.Data.eulerangles;

