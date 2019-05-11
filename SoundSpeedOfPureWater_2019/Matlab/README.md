# Requirements:
* MATLAB R2016a or higher, to support the demo scripts.  In versions R2014b or higher, live scripts will open, but with the code only (no formatted text or in-line output).
* Curve Fitting Toolbox
  * [Available for a 30-day trial if not installed](https://www.mathworks.com/campaigns/products/trials.html?prodcode=CF)
  * Since the LBFs use tensor spline basis functions, two functions from the toolbox are required:
    * **fnval**: evaluates a spline either from scattered data or gridded data
    * **fnder**: returns a specified derivative of a spline representation
* Custom MATLAB functions, available as a snapshot version from this archive, or the latest version in the [jmichaelb/LocalBasisFunction Github repository](https://github.com/jmichaelb/LocalBasisFunction/tree/master/Matlab)
  * **fnGval**: evaluates thermodynamic properties using the LBF for Gibbs Energy
  * **therm_surf**: creates 3D thermmdynamic surfaces as a function of P and T
  * **IAPWS95**: returns properties based on the IAPWS-95 representation
* In the MATLAB path, a data file WaterEOS.mat containing both a collection of relevant data and the LBF representation.


