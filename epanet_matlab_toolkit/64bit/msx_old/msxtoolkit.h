/******************************************************************************
**  MODULE:        MSXTOOLKIT.H
**  PROJECT:       EPANET-MSX
**  DESCRIPTION:   C/C++ header file for MSX API toolkit.
**  COPYRIGHT:     Copyright (C) 2007 Feng Shang, Lewis Rossman, and James Uber.
**                 All Rights Reserved. See license information in LICENSE.TXT.
**  AUTHORS:       K. Arrowood, Xylem intern
**  VERSION:       1.1 
**  LAST UPDATE:   Refer to git history
*******************************************************************************/
 
#ifndef ENUMSOPEN
#include "msxenums.h"
#endif

// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

// --- define DLLEXPORT

#undef DLLEXPORT
#ifdef WINDOWS
  #ifdef __cplusplus
  #define DLLEXPORT extern "C" __declspec(dllexport) __stdcall
  #else
  #define DLLEXPORT __declspec(dllexport) __stdcall
  #endif
#else
  #ifdef __cplusplus
  #define DLLEXPORT extern "C"
  #else
  #define DLLEXPORT
  #endif
#endif





//Project Functions
int DLLEXPORT MSXopen(void);
int DLLEXPORT MSXclose(void);
int DLLEXPORT MSXinit(void);
int DLLEXPORT MSXprintQuality(int type, char *id, char *species, char *fname);


//Network building functions
int DLLEXPORT MSXaddNode(char *id);
int DLLEXPORT MSXaddTank(char *id, double initialVolume, int mixModel, double volumeMix);
int DLLEXPORT MSXaddReservoir(char *id, double initialVolume, int mixModel, double volumeMix);
int DLLEXPORT MSXaddLink(char *id, char *startNode, char *endNode, double length, double diameter, double roughness);

//Species/Chemistry option functions
int DLLEXPORT MSXaddOption(int optionType, char * value);
int DLLEXPORT MSXaddSpecies(char *id, int type, int units, double aTol, double rTol);
int DLLEXPORT MSXaddCoefficeint(int type, char *id, double value);
int DLLEXPORT MSXaddTerm(char *id, char *equation);
int DLLEXPORT MSXaddExpression(int classType, int expressionType, char *species, char *equation);
int DLLEXPORT MSXaddSource(int sourceType, char *nodeId, char *speciesId, double strength, char *timePattern);
int DLLEXPORT MSXaddQuality(char *type, char *speciesId, double value, char *id);
int DLLEXPORT MSXaddParameter(char *type, char *paramId, double value, char *id);
int DLLEXPORT MSXsetReport(char *reportType, char *id, int precision);

//Hydraulic Functions
int DLLEXPORT MSXsetHydraulics(float *demands, float *heads, float *flows);

int DLLEXPORT MSXsetSize(int type, int size);


// Below is from the legacy epanetmsx.h

int  DLLEXPORT MSXgetindex(int type, char *id, int *index);
int  DLLEXPORT MSXgetIDlen(int type, int index, int *len);
int  DLLEXPORT MSXgetID(int type, int index, char *result, int len);
int  DLLEXPORT MSXgetcount(int type, int *count);
int  DLLEXPORT MSXgetspecies(int index, int *type, char *units, double *aTol,
               double *rTol);
int  DLLEXPORT MSXgetconstant(int index, double *value);
int  DLLEXPORT MSXgetparameter(int type, int index, int param, double *value);
int  DLLEXPORT MSXgetsource(int node, int species, int *type, double *level,
               int *pat);
int  DLLEXPORT MSXgetpatternlen(int pat, int *len);
int  DLLEXPORT MSXgetpatternvalue(int pat, int period, double *value);
int  DLLEXPORT MSXgetinitqual(int type, int index, int species, double *value);
int  DLLEXPORT MSXgetQualityByIndex(int type, int index, int species, double *value);
int  DLLEXPORT MSXgetQualityByID(int type, char *id, char *species, double *value);
int  DLLEXPORT MSXsetconstant(int index, double value);
int  DLLEXPORT MSXsetparameter(int type, int index, int param, double value);
int  DLLEXPORT MSXsetinitqual(int type, int index, int species, double value);
int  DLLEXPORT MSXsetsource(int node, int species, int type, double level,
               int pat);
int  DLLEXPORT MSXsetpatternvalue(int pat, int period, double value);
int  DLLEXPORT MSXaddpattern(char *id);
int  DLLEXPORT MSXsetpattern(int pat, double mult[], int len);

int DLLEXPORT MSXstep(long *t, long *tleft);
int DLLEXPORT MSXgeterror(int code, char *msg, int len);

//Simulation Options
int DLLEXPORT MSXsetFlowFlag(int flag);
int DLLEXPORT MSXsetTimeParameter(int type, long value);