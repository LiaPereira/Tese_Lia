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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

//espetro
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
//fim espetro

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),spectrum("off")
{
  G4int n_particle = 1;
  //fpParticleGun  = new G4ParticleGun(n_particle);
  fpParticleGun  = new G4GeneralParticleSource();
  // default particle kinematic

  //  G4ParticleDefinition* particle
  //  = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  //  fpParticleGun->SetParticleDefinition(particle);
  //  fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //  //fpParticleGun->SetParticleEnergy(0.1*MeV);
  //  fpParticleGun->SetParticlePosition(G4ThreeVector(0.*nm,0.*nm,0.*nm));

  // ELE ENTRA AQUI SÓ UMA VEZ (NO INICIO DO RUN) E TEM A ENERGIA, POSIÇÃO, ETC COMO DEFAULT

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fpParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //espetro
   if (spectrum =="on")
    {
      G4cout<<"hey"<<G4endl;
      // isto em baixo é uma sequencial runaction (????) - G4RunManager
      // Há a opção de MT master runaction com o G4MTRunManager
    runAction = static_cast<const RunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()); 

	  G4DataVector* energies =  runAction->GetEnergies();
	  G4DataVector* data =  runAction->GetData();
    //G4cout<<"Energias " <<*energies<<G4endl;
	 
	  G4double sum = runAction->GetDataSum();
    //G4cout<<"A soma é: " <<sum<<G4endl;
	  G4double partSum = 0;
	  G4int j = 0;
	  G4double random= sum*G4UniformRand();
    //G4cout<<"Random: " <<random<<G4endl;
	  while (partSum<random)
	    {
	      partSum += (*data)[j];
	      j++;    
        
	    }
    G4cout<<"energia "<<(*energies)[j-1]<<G4endl;
	 
	  //fpParticleGun->SetParticleEnergy((*energies)[j-1]); //Eu acrescentei o -1 !!!!!!!!!!!!!!
  //fim espetro
	
    }

  G4cout<<"Energia da particula: " <<fpParticleGun->GetParticleEnergy()<<G4endl; 
  G4cout<<"Posição da particula: " <<fpParticleGun->GetParticlePosition()<<G4endl; 
  G4cout<<"Direção da particula: " <<fpParticleGun->GetParticleMomentumDirection()<<G4endl; 
  G4cout<<"Definição da particula: " <<fpParticleGun->GetParticleDefinition()->GetParticleName()<<G4endl; 
  

  
  fpParticleGun->GeneratePrimaryVertex(anEvent);
}
