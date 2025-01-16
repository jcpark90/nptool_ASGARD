#ifndef TWAS3ABISPECTRA_H
#define TWAS3ABISPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta                 address: a.matta@surrey.ac.uk   *
 *                                                                           *
 * Creation Date  : June 2015                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for WAS3ABi                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TWAS3ABiData.h"
#include "TWAS3ABiPhysics.h"

// ForwardDeclaration
class TWAS3ABiPhysics ;

class TWAS3ABiSpectra:public VSpectra {
  public:
    // constructor and destructor
    TWAS3ABiSpectra();
    TWAS3ABiSpectra(unsigned int NumberOfDetector);
    ~TWAS3ABiSpectra();

  private:
    // Initialization methods
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  public:
    // Filling methods
    void FillRawSpectra(TWAS3ABiData*);
    void FillPreTreatedSpectra(TWAS3ABiData*);
    void FillPhysicsSpectra(TWAS3ABiPhysics*);

  private: // Information on WAS3ABI
    unsigned int fNumberOfDetector;
    unsigned int fStripFront;
    unsigned int fStripBack;
    unsigned int fPad;
};

#endif
