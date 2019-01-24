/******************************************************************************
**  MODULE:        EPANETMSX.H
**  PROJECT:       EPANET-MSX
**  DESCRIPTION:   C/C++ header file for EPANET Multi-Species Extension Toolkit
**  COPYRIGHT:     Copyright (C) 2007 Feng Shang, Lewis Rossman, and James Uber.
**                 All Rights Reserved. See license information in LICENSE.TXT.
**  AUTHORS:       L. Rossman, US EPA - NRMRL
**                 F. Shang, University of Cincinnati
**                 J. Uber, University of Cincinnati
**  VERSION:       1.1.00 
**  LAST UPDATE:   7/31/07
*******************************************************************************/

#ifndef EPANETMSX_H
#define EPANETMSX_H

// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

// --- define DLLEXPORT
#ifndef DLLEXPORT
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
#endif
// --- define MSX constants

#define MSX_NODE      0
#define MSX_LINK      1
#define MSX_TANK      2
#define MSX_SPECIES   3
#define MSX_TERM      4
#define MSX_PARAMETER 5
#define MSX_CONSTANT  6
#define MSX_PATTERN   7

#define MSX_BULK      0
#define MSX_WALL      1

#define MSX_NOSOURCE  -1
#define MSX_CONCEN     0
#define MSX_MASS       1
#define MSX_SETPOINT   2
#define MSX_FLOWPACED  3

// --- declare MSX functions

int  DLLEXPORT MSXopen(char *fname);
int  DLLEXPORT MSXsolveH(void);
int  DLLEXPORT MSXusehydfile(char *fname);
int  DLLEXPORT MSXsolveQ(void);
int  DLLEXPORT MSXinit(int saveFlag);
int  DLLEXPORT MSXstep(long *t, long *tleft);
int  DLLEXPORT MSXsaveoutfile(char *fname);
int  DLLEXPORT MSXsavemsxfile(char *fname);
int  DLLEXPORT MSXreport(void);
int  DLLEXPORT MSXclose(void);

int  DLLEXPORT MSXgetindex(int type, char *id, int *index);
int  DLLEXPORT MSXgetIDlen(int type, int index, int *len);
int  DLLEXPORT MSXgetID(int type, int index, char *id, int len);
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
int  DLLEXPORT MSXgetqual(int type, int index, int species, double *value);
int  DLLEXPORT MSXgeterror(int code, char *msg, int len);

int  DLLEXPORT MSXsetconstant(int index, double value);
int  DLLEXPORT MSXsetparameter(int type, int index, int param, double value);
int  DLLEXPORT MSXsetinitqual(int type, int index, int species, double value);
int  DLLEXPORT MSXsetsource(int node, int species, int type, double level,
               int pat);
int  DLLEXPORT MSXsetpatternvalue(int pat, int period, double value);
int  DLLEXPORT MSXsetpattern(int pat, double mult[], int len);
int  DLLEXPORT MSXaddpattern(char *id);

#endif
