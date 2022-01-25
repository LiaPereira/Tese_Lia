
#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
  :Action(Gun)
{ 

spectrum = new G4UIcmdWithAString("/gun/spectrum",this);
  spectrum->SetGuidance("Shoot the incident particle with a certain energy spectrum.");
  spectrum->SetGuidance("  Choice : on(default), off");
  spectrum->SetParameterName("choice",true);
  spectrum->SetDefaultValue("on");
  spectrum->SetCandidates("on off");
  spectrum->AvailableForStates(G4State_PreInit,G4State_Idle);

}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{

delete spectrum;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 

if( command == spectrum )
{ Action->SetSpectrum(newValue);} 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
