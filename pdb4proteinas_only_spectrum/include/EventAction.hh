//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

#include <map>

class EventActionMessenger;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void AddEdepEvent(G4double edep)
  {
    fTotalEnergyDeposit += edep;
  };
  G4double GetEdepEvent()
  {
    return fTotalEnergyDeposit;
  };

  // ---------------- EU ---------------------------------------------
  void AddEdepToResidue(G4int numMolec,G4int numResi,G4double edep)
  {
    if(numMolec==1)
    {
      fEdepMolec1[numResi]+=edep;
    }
    else{
      fEdepMolec2[numResi]+=edep;
    }
  }
  // ----------------------------------------------------------------
  void SetPrintModulo(G4int val)
  {
    fPrintModulo = val;
  };


private:
  // total energy deposit per event
  G4double fTotalEnergyDeposit;

  // ------------------ EU ------------------
  std::map<G4int,G4double> fEdepMolec1;
  std::map<G4int,G4double> fEdepMolec2;

  // map: (G4int : nucleotide ID, G4double : energy deposit)
  // -----------------------------------------

  G4int                     fPrintModulo;
  EventActionMessenger*     fpEventMessenger;

  // Compute Strand breaks from energy deposits in the molecules

  //  ------------- EU ------------------------
  void ComputeMoleculeBreaks(G4int*);
  // ------------------------------------------
};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
