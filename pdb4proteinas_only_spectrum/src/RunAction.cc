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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "RunInitObserver.hh"

#include "Analysis.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"

#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

//espetro
#include "G4DataVector.hh"
#include "G4SystemOfUnits.hh"
//fim espetro

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction() : G4UserRunAction()
{

  i = 0;

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  i++;
  RunInitManager::Instance()->Initialize();

  //espetro
  energies = new G4DataVector;
  data = new G4DataVector;
  
  G4cout<<"----- Vai chamar a função ReadData -----"<<G4endl;
  //ReadData(keV,"M_flare");
  //fim espetro

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFirstHistoId(1);

  // Creating histograms
  //
  analysisManager->CreateH1("1",
                            "Energy deposit in the target (eV)",
                            1000.,0.,1000.);
  analysisManager->CreateH1("2",
                            "Number of SSB",
                            10,0.,10.);
  analysisManager->CreateH1("3",
                            "Number of DSB",
                            10,0.,10.);

  // Open an output file
  //

  // vários outputs

  const PrimaryGeneratorAction *primaryGeneratorAction = static_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  
  const G4GeneralParticleSource *fpParticleGun = primaryGeneratorAction->GetParticleGun();
  G4String Energia = G4UIcommand::ConvertToString(fpParticleGun->GetParticleEnergy()/CLHEP::keV) + "keV" ;

  //G4cout<<"Energia do feixe "<<Energia<<G4endl; //está mal
  
  G4String i_ = G4UIcommand::ConvertToString(i);

  G4String fileName = "pdb4dna_output_"+i_;
  analysisManager->OpenFile(fileName);
  G4String extension = analysisManager->GetFileType();
  fileName = fileName + "." + extension;

  G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // save histograms 
  //
  analysisManager->Write();

  analysisManager->CloseFile();

  delete G4AnalysisManager::Instance();
}

//espetro
void RunAction::ReadData(G4double unitE, G4String fileName)
{
  std::ostringstream ost;
  
  ost << fileName <<".dat";
  
  G4String name = ost.str();
  char* path;

  G4cout<<"O nome do ficheiro é "<< name<<G4endl;
  
  if (!(std::getenv("XRAYDATA"))) { 
    // Ele entra aqui
    path = std::getenv("PWD");    
  }
  
  else {    
    path = std::getenv("XRAYDATA");
  }
  
  
  G4String pathString(path);
  name = pathString + "/" + name;
  G4cout<<"O path do ficheiro é "<< name<<G4endl;
  
  std::ifstream file(name);
  std::filebuf* lsdp = file.rdbuf();

  
  if (! (lsdp->is_open()) )
    {
      // Ele nao entra aqui (so deve entrar se nao encontrar o ficheiro)
      G4ExceptionDescription execp;
      execp <<  "XrayFluoRunAction - data file: " + name + " not found";
      G4Exception("XrayFluoRunAction::ReadData()","example-xray_fluorescence04",
	  FatalException, execp);
    }
  G4double a = 0;
  G4int k = 1;
  G4cout<<"Vamos percorrer as linhas e as colunas do ficheiro:" <<G4endl;
  // O "do" só acaba no while lá em baixo
  do
    {
      file >> a;
      //G4cout<<a<<G4endl;
      G4int nColumns = 2;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      if (a == 200 || a == 200)
	{
	  
	}
      else
	{
	  if (k%nColumns != 0)
	    {	
        // Entra neste if quando o "a" corresponde a uma energia (1a coluna)
        // As energias ficam guardadas
	      G4double e = a * unitE;
	      energies->push_back(e);
	      
	      k++;
	      
	    }
	  else if (k%nColumns == 0)
	    {
        // Entra neste if quando o "a" corresponde a uma "fluencia/quantidade" (2a coluna)
        // Os dados correspondentes a cada energia ficam guardados
	      G4double value = a;
	      data->push_back(value);
	      
	      k = 1;
	    }
	}
      
    } while (a != 20); // end of file
  
  file.close();
  G4cout << " done" << G4endl;
}

G4DataVector* RunAction::GetEnergies() const
{
  return energies;
}

G4DataVector* RunAction::GetData() const
{
  return data;
}

G4double RunAction::GetDataSum() const
{
 
  G4double sum = 0;
  for (size_t i = 0; i < data->size(); i++)
    {
      sum+=(*data)[i];
    }
  return sum;
}
//fim espetro