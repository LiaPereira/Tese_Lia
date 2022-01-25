#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
  ~PrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);

  private: 

  PrimaryGeneratorAction* Action;

  G4UIcmdWithAString*  spectrum;

  };



#endif