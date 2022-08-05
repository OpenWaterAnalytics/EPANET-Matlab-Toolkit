/******************************************************************************
**  MODULE:        LEGACYTOOLKIT.H
**  PROJECT:       EPANET-MSX
**  DESCRIPTION:   C/C++ header file for EPANET Multi-Species Extension Toolkit
**  COPYRIGHT:     Copyright (C) 2007 Feng Shang, Lewis Rossman, and James Uber.
**                 All Rights Reserved. See license information in LICENSE.TXT.
**  AUTHORS:       L. Rossman, US EPA - NRMRL
**                 F. Shang, University of Cincinnati
**                 J. Uber, University of Cincinnati
**                 K. Arrowood, Xylem intern
**  VERSION:       1.1 
**  LAST UPDATE:   Refer to git history
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

// --- declare MSX functions

int  DLLEXPORT Legacyopen(MSXproject *MSX, char*argv[]);
int  DLLEXPORT MSXsolveH(MSXproject MSX);
int  DLLEXPORT MSXusehydfile(MSXproject MSX);
int  DLLEXPORT MSXsolveQ(MSXproject MSX);
int  DLLEXPORT Legacyinit(MSXproject MSX, int saveFlag);
int  DLLEXPORT MSXsaveoutfile(MSXproject MSX, char *fname);
int  DLLEXPORT MSXsavemsxfile(MSXproject MSX, char *fname);
int  DLLEXPORT MSXreport(MSXproject MSX, char *fname);
int  DLLEXPORT Legacyclose(MSXproject MSX);

int  DLLEXPORT MSXsaveResults(MSXproject MSX);
int  DLLEXPORT MSXsaveFinalResults(MSXproject MSX);

void DLLEXPORT MSXrunLegacy(int argc, char *argv[]);


#endif