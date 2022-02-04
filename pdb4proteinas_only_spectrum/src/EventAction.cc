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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "Analysis.hh"
#include "EventActionMessenger.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>

#include "G4RunManager.hh"
#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction():G4UserEventAction()
{
  //default parameter values
  //
  fPrintModulo=100;
  //fThresEdepForSSB=8.22*eV;
  //fThresDistForDSB=10;
  fTotalEnergyDeposit=0;

  //create commands
  //
  fpEventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete fpEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oo ooOO0OOooo......

void EventAction::BeginOfEventAction( const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();


  // Initialization of parameters
  //
  fTotalEnergyDeposit=0.;


  //  ------------------ EU -------------------------
  fEdepMolec1.clear();
  fEdepMolec2.clear();
  list_resi.clear();
  list_elem.clear();
  list_numelem.clear();
  // ----------------------------------
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction( const G4Event*)
{
  // At the end of an event, compute the number of strand breaks
  //
  
  //  ----------------- EU -------------------
  G4int mb[1] = {0};
  ComputeMoleculeBreaks(mb);
  // -----------------------------------------

  // Fill histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  std::ofstream file;
  
  if ( fTotalEnergyDeposit>0. )
  {
    analysisManager->FillH1(1,fTotalEnergyDeposit);

    // ------------ FICHEIRO TXT --------------------------

    const RunAction *runAction = static_cast<const RunAction*>(G4RunManager::GetRunManager()->GetUserRunAction());
    G4String run = G4UIcommand::ConvertToString(runAction -> Run);
    std::string txt_file = "file_run_"+run;
    G4cout<<"hello"<<G4endl;

    file.open(txt_file, std::ios_base::app | std::ios_base::in);
    if (file.is_open())
    {
      
      if (!list_resi.empty())
      {       

        //loop para anotarmos todos os resíduos (e elementos correspondentes) de cada evento

        //file << fTotalEnergyDeposit << "\t" << mb[0] << "\t"<<list_resi[0]<< "\t"<<list_elem[0][0]<< "\t"<<list_numelem[0]<<G4endl;
          // G4int i;
          // G4int k;
          // for(i = 1;i<list_resi.size();i++)
          // {   file << "----" << "\t" << "--" << "\t"<<list_resi[i]<< "\t"<<list_elem[1][0]<< "\t"<<list_numelem[0]<<G4endl;
          //     for (k=1;k<list_elem[i].size();k++)
          //     {
          //       file << "----" << "\t" << "--" << "\t"<<"--"<< "\t"<<list_elem[i][k]<< "\t"<<list_numelem[0]<<G4endl;
          //     }
          // }

        // A AUTENTICA TRAPALHADA:

        G4int i;
        G4int k;
        std::string residuos = "[";
        std::string elementos = "[";
        std::string num_elementos = "[";
        for (i=0;i<list_resi.size();i++)
        {
  
          residuos.append(std::to_string(list_resi[i]));

          if (i!=list_resi.size()-1)
          {residuos.append(",");}
          

          elementos.append("(");
          num_elementos.append("(");

          for (k=0;k<list_elem[i].size();k++)
          {
             elementos.append(list_elem[i][k]);
             num_elementos.append(std::to_string(list_numelem[i][k]));
          
             if (k!=list_elem[i].size()-1)
             {
               elementos.append(",");
               num_elementos.append(",");
             }
          }

          elementos.append(")");
          num_elementos.append(")");
          if (i!=list_elem.size()-1)
          {
            elementos.append(",");
            num_elementos.append(",");
          }

        }

        residuos.append("]");
        elementos.append("]");
        num_elementos.append("]");

        file << fTotalEnergyDeposit << "\t" << mb[0] << "\t"<<residuos<< "\t"<<elementos<< "\t"<<num_elementos<<G4endl;  

      //AQUI ACABA A TRAPALHADA
        

      }
      else
      {
      file << fTotalEnergyDeposit << "\t" << mb[0] <<G4endl;
      }
    }
    file.close();
    // -------------------------------------------------------
      
  }

  //  ------------- EU -------------------
  if (mb[0]>0)
  {
  analysisManager->FillH1(2, mb[0]);
  // -------------------------------------
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//  -------- EU ----------------------------------------------
void EventAction::ComputeMoleculeBreaks(G4int* mb)
{
  // sb quantities
  //
  G4int mb_=0;

  // residue id and energy deposit
  G4int resi;
  G4double edep;
  
  while ( !fEdepMolec1.empty() )
  {
    resi = fEdepMolec1.begin()->first;
    edep = fEdepMolec1.begin()->second;
    fEdepMolec1.erase( fEdepMolec1.begin() );
    //G4cout<<"edep "<<edep<<G4endl;

    if ( edep >= 0/eV ) // Estou confusa com as unidades do edep. 
    {
      mb_++;
      //G4cout<<"É maior que 0.09"<<G4endl;
    }
  
  }
  mb[0]=mb_;
}
// -----------------------------------------------------------------------
