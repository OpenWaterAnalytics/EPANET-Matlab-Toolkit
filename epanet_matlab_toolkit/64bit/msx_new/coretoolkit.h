/******************************************************************************
**  MODULE:        CORETOOLKIT.H
**  PROJECT:       EPANET-MSX
**  DESCRIPTION:   C/C++ header file for new MSX core API toolkit.
**  COPYRIGHT:     Copyright (C) 2007 Feng Shang, Lewis Rossman, and James Uber.
**                 All Rights Reserved. See license information in LICENSE.TXT.
**  AUTHORS:       L. Rossman, US EPA - NRMRL
**                 F. Shang, University of Cincinnati
**                 J. Uber, University of Cincinnati
**                 K. Arrowood, Xylem intern
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

// Opaque Pointer
typedef struct Project *MSXproject;

//Project Functions
int DLLEXPORT MSX_open(MSXproject *MSX);
int DLLEXPORT MSX_close(MSXproject MSX);
int DLLEXPORT MSX_init(MSXproject MSX);
int DLLEXPORT MSX_printQuality(MSXproject MSX, int type, char *id, char *species, char *fname);


//Network building functions
int DLLEXPORT MSX_addNode(MSXproject MSX, char *id);
int DLLEXPORT MSX_addTank(MSXproject MSX, char *id, double initialVolume, int mixModel, double volumeMix);
int DLLEXPORT MSX_addReservoir(MSXproject MSX, char *id, double initialVolume, int mixModel, double volumeMix);
int DLLEXPORT MSX_addLink(MSXproject MSX, char *id, char *startNode, char *endNode, double length, double diameter, double roughness);

//Species/Chemistry option functions
int DLLEXPORT MSX_addOption(MSXproject MSX, int optionType, char * value);
int DLLEXPORT MSX_addSpecies(MSXproject MSX, char *id, int type, int units, double aTol, double rTol);
int DLLEXPORT MSX_addCoefficeint(MSXproject MSX, int type, char *id, double value);
int DLLEXPORT MSX_addTerm(MSXproject MSX, char *id, char *equation);
int DLLEXPORT MSX_addExpression(MSXproject MSX, int classType, int expressionType, char *species, char *equation);
int DLLEXPORT MSX_addSource(MSXproject MSX, int sourceType, char *nodeId, char *speciesId, double strength, char *timePattern);
int DLLEXPORT MSX_addQuality(MSXproject MSX, char *type, char *speciesId, double value, char *id);
int DLLEXPORT MSX_addParameter(MSXproject MSX, char *type, char *paramId, double value, char *id);
int DLLEXPORT MSX_setReport(MSXproject MSX, char *reportType, char *id, int precision);

//Hydraulic Functions
int DLLEXPORT MSX_setHydraulics(MSXproject MSX, float *demands, float *heads, float *flows);

int DLLEXPORT MSX_setSize(MSXproject MSX, int type, int size);


// Below is from the legacy epanetmsx.h

int  DLLEXPORT MSX_getindex(MSXproject MSX, int type, char *id, int *index);
int  DLLEXPORT MSX_getIDlen(MSXproject MSX, int type, int index, int *len);
int  DLLEXPORT MSX_getID(MSXproject MSX, int type, int index, char *id, int len);
int  DLLEXPORT MSX_getcount(MSXproject MSX, int type, int *count);
int  DLLEXPORT MSX_getspecies(MSXproject MSX, int index, int *type, char *units, double *aTol,
               double *rTol);
int  DLLEXPORT MSX_getconstant(MSXproject MSX, int index, double *value);
int  DLLEXPORT MSX_getparameter(MSXproject MSX, int type, int index, int param, double *value);
int  DLLEXPORT MSX_getsource(MSXproject MSX, int node, int species, int *type, double *level,
               int *pat);
int  DLLEXPORT MSX_getpatternlen(MSXproject MSX, int pat, int *len);
int  DLLEXPORT MSX_getpatternvalue(MSXproject MSX, int pat, int period, double *value);
int  DLLEXPORT MSX_getinitqual(MSXproject MSX, int type, int index, int species, double *value);
int  DLLEXPORT MSX_getQualityByIndex(MSXproject MSX, int type, int index, int species, double *value);
int  DLLEXPORT MSX_getQualityByID(MSXproject MSX, int type, char *id, char *species, double *value);
int  DLLEXPORT MSX_setconstant(MSXproject MSX, int index, double value);
int  DLLEXPORT MSX_setparameter(MSXproject MSX, int type, int index, int param, double value);
int  DLLEXPORT MSX_setinitqual(MSXproject MSX, int type, int index, int species, double value);
int  DLLEXPORT MSX_setsource(MSXproject MSX, int node, int species, int type, double level,
               int pat);
int  DLLEXPORT MSX_setpatternvalue(MSXproject MSX, int pat, int period, double value);
int  DLLEXPORT MSX_addpattern(MSXproject MSX, char *id);
int  DLLEXPORT MSX_setpattern(MSXproject MSX, int pat, double mult[], int len);

int DLLEXPORT MSX_step(MSXproject MSX, long *t, long *tleft);

//Simulation Options
int DLLEXPORT MSX_setFlowFlag(MSXproject MSX, int flag);
int DLLEXPORT MSX_setTimeParameter(MSXproject MSX, int type, long value);