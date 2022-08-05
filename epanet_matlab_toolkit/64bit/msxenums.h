/************************************************************************
**  MODULE:        MSXENUMS.H
**  PROJECT:       EPANET-MSX
**  DESCRIPTION:   Global enums for MSX toolkit.
**  COPYRIGHT:     Copyright (C) 2007 Feng Shang, Lewis Rossman, and James Uber.
**                 All Rights Reserved. See license information in LICENSE.TXT.
**  AUTHORS:       L. Rossman, US EPA - NRMRL
**                 F. Shang, University of Cincinnati
**                 J. Uber, University of Cincinnati
**                 K. Arrowood, Xylem intern
**  VERSION:       1.1.00
**  LAST UPDATE:   Refer to git history
***********************************************************************/

#define ENUMSOPEN 1

//-----------------------------------------------------------------------------
//  Enumerated Types
//-----------------------------------------------------------------------------
 enum ObjectType                       // Object types                         //1.1.00
                {NODE,
                 LINK,
                 TANK,
                 SPECIES,
                 TERM,
                 PARAMETER,
                 CONSTANT,
                 PATTERN,
                 MAX_OBJECTS,
                 };

 enum UnitSystemType                   // Unit system:
                 {US,                  //   US
                  SI};                 //   SI (metric)

 enum FlowUnitsType                    // Flow units:
                 {CFS,                 //   cubic feet per second
                  GPM,                 //   gallons per minute
                  MGD,                 //   million gallons per day
                  IMGD,                //   imperial million gal. per day
                  AFD,                 //   acre-feet per day
                  LPS,                 //   liters per second
                  LPM,                 //   liters per minute
                  MLD,                 //   megaliters per day
                  CMH,                 //   cubic meters per hour
                  CMD};                //   cubic meters per day

 enum MixType                          // Tank mixing regimes
                 {MIX1,                //   1-compartment model
                  MIX2,                //   2-compartment model
                  FIFO,                //   First in, first out model
                  LIFO};               //   Last in, first out model

 enum SpeciesType                      // Types of water quality species
                {BULK,                 //   bulk flow species
                 WALL};                //   pipe wall attached species

 enum ExpressionType                   // Types of math expressions
                {NO_EXPR,              //   no expression
                 RATE,                 //   reaction rate
                 FORMULA,              //   simple formula
                 EQUIL};               //   equilibrium expression

 enum SolverType                       // ODE solver options
                 {EUL,                 //   Euler
                  RK5,                 //   5th order Runge-Kutta
                  ROS2};               //   2nd order Rosenbrock

 enum CouplingType                     // Degree of coupling for solving DAE's
                 {NO_COUPLING,         //   no coupling between alg. & diff. eqns.
                  FULL_COUPLING};      //   full coupling between alg. &diff. eqns.

 enum MassUnitsType                    // Concentration mass units
                 {MG,                  //   milligram
                  UG,                  //   microgram
                  MOLE,                //   mole
                  MMOLE};              //   millimole

 enum AreaUnitsType                    // Pipe surface area units
                 {FT2,                 //   square feet
                  M2,                  //   square meters
                  CM2};                //   square centimeters

 enum RateUnitsType                    // Reaction rate time units
                 {SECONDS,             //   seconds
                  MINUTES,             //   minutes
                  HOURS,               //   hours
                  DAYS};               //   days

 enum UnitsType                        // Measurement unit types
                 {LENGTH_UNITS,        //   length
                  DIAM_UNITS,          //   pipe diameter
                  AREA_UNITS,          //   surface area
                  VOL_UNITS,           //   volume
                  FLOW_UNITS,          //   flow
                  CONC_UNITS,          //   concentration volume
                  RATE_UNITS,          //   reaction rate time units
                  MAX_UNIT_TYPES};

 enum HydVarType                       // Hydraulic variables
                 {DIAMETER = 1,        //   link diameter
                  FLOW,                //   link flow rate
                  VELOCITY,            //   link flow velocity
                  REYNOLDS,            //   Reynolds number
                  SHEAR,               //   link shear velocity
                  FRICTION,            //   friction factor
                  AREAVOL,             //   area/volume
                  ROUGHNESS,           //   roughness                          /*Feng Shang 01/29/2008*/
                  MAX_HYD_VARS};

 enum TstatType                        // Time series statistics
                 {SERIES,              //   full time series
                  AVGERAGE,            //   time-averages
                  MINIMUM,             //   minimum values
                  MAXIMUM,             //   maximum values
                  RANGE};              //   max - min values

 enum OptionType                       // Analysis options
                 {AREA_UNITS_OPTION,
                  RATE_UNITS_OPTION,
                  SOLVER_OPTION,
                  COUPLING_OPTION,
                  TIMESTEP_OPTION,
                  RTOL_OPTION,
                  ATOL_OPTION,
                  COMPILER_OPTION};                                            //1.1.00

 enum CompilerType                     // C compiler type                      //1.1.00
                 {NO_COMPILER,
                  VC,                  // MS Visual C compiler
                  GC};                 // Gnu C compiler

 enum FileModeType                     // File modes
                 {SCRATCH_FILE,
                  SAVED_FILE,
                  USED_FILE};

 enum SectionType                      // Input data file sections
                 {s_TITLE,
                  s_SPECIES,
                  s_COEFF,
                  s_TERM,
                  s_PIPE,
                  s_TANK,
                  s_SOURCE,
                  s_QUALITY,
                  s_PARAMETER,
                  s_PATTERN,
                  s_OPTION,
                  s_REPORT};

 enum InpErrorCodes {                   // Error codes (401 - 409)
        INP_ERR_FIRST        = 400,
        ERR_LINE_LENGTH,
        ERR_ITEMS, 
        ERR_KEYWORD,
        ERR_NUMBER,
        ERR_NAME,
        ERR_RESERVED_NAME,
        ERR_DUP_NAME,
        ERR_DUP_EXPR, 
        ERR_MATH_EXPR,
        INP_ERR_LAST};
 
 
 enum ErrorCodeType                    // Error codes (501-524)
          {ERR_FIRST = 500,
           ERR_MEMORY,                 // 501
           ERR_NO_EPANET_FILE,         // 502
           ERR_OPEN_MSX_FILE,          // 503
           ERR_OPEN_HYD_FILE,          // 504
           ERR_READ_HYD_FILE,          // 505
           ERR_MSX_INPUT,              // 506
           ERR_NUM_PIPE_EXPR,          // 507
           ERR_NUM_TANK_EXPR,          // 508
           ERR_INTEGRATOR_OPEN,        // 509
           ERR_NEWTON_OPEN,            // 510
           ERR_OPEN_OUT_FILE,          // 511
           ERR_IO_OUT_FILE,            // 512
           ERR_INTEGRATOR,             // 513
           ERR_NEWTON,                 // 514
           ERR_INVALID_OBJECT_TYPE,    // 515
           ERR_INVALID_OBJECT_INDEX,   // 516
           ERR_UNDEFINED_OBJECT_ID,    // 517
           ERR_INVALID_OBJECT_PARAMS,  // 518
           ERR_MSX_NOT_OPENED,         // 519
           ERR_MSX_OPENED,             // 520
           ERR_OPEN_RPT_FILE,          // 521                                  //(LR-11/20/07, to fix bug 08)
           ERR_COMPILE_FAILED,         // 522                                  //1.1.00
           ERR_COMPILED_LOAD,          // 523                                  //1.1.00
		       ERR_ILLEGAL_MATH,           // 524                                  //1.1.00
           ERR_HYD,
           ERR_INIT,
           ERR_MAX};

/// Time parameters (From EPANET)
/**
These time-related options are used with gettimeparam and settimeparam.
All times are expressed in seconds The parameters marked as read only are
computed values that can only be retrieved.
*/
typedef enum {
  DURATION,  //!< Total simulation duration
  HYDSTEP,  //!< Hydraulic time step
  QUALSTEP,  //!< Water quality time step
  PATTERNSTEP,  //!< Time pattern period
  PATTERNSTART,  //!< Time when time patterns begin
  REPORTSTEP,  //!< Reporting time step
  REPORTSTART,  //!< Time when reporting starts
  STATISTIC,  //!< Reporting statistic code
  HTIME, //!< Elapsed time of current hydraulic solution (read only)
  QTIME, //!< Elapsed time of current quality solution (read only)
} TimeParameter;

/// Size Limts
/**
Limits on the size of character arrays used to store ID names
and text messages.
*/
typedef enum {
  MAXID   = 31,     //!< Max. # characters in ID name
  MAXMSG  = 255     //!< Max. # characters in message text
} SizeLimits;

typedef enum {  // Type of source quality input
  NOSOURCE  = -1,
  CONCEN    = 0,  //    inflow concentration
  MASS      = 1,  //    mass inflow booster
  SETPOINT  = 2,  //    setpoint booster
  FLOWPACED = 3   //    flow paced booster
} SourceTypes;
     