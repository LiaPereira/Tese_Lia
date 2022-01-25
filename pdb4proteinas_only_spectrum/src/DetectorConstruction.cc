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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fpMessenger(0),
  fCheckOverlaps(false)
{
  //Select a PDB file name by default
  //otherwise modified by the LoadPDBfile messenger
  //
  fPdbFileName=G4String("1iyt.pdb");
  fPdbFileStatus=0;
  fChosenOption=11;

  fpDefaultMaterial=0;
  fpWaterMaterial=0;
  fpMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  fChosenOption=11;//Draw atomic with bounding box (visualisation by default)
  fPdbFileStatus=0;//There is no PDB file loaded at this stage

  //Define materials and geometry
  G4VPhysicalVolume* worldPV;
  worldPV=DefineVolumes(fPdbFileName,fChosenOption);
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
  //[G4_WATER]
  //
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial("G4_WATER", fromIsotopes);
  fpWaterMaterial = G4Material::GetMaterial("G4_WATER");

  //[Vacuum]
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                 kStateGas, 2.73*kelvin, 3.e-18*pascal);
  fpDefaultMaterial = G4Material::GetMaterial("Galactic");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CheckMaterials()
{
  if ( !fpDefaultMaterial || !fpWaterMaterial )
  {
    G4cerr << "Cannot retrieve materials already defined. " << G4endl;
    G4cerr << "Exiting application " << G4endl;
    exit(1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructWorld()
{
  // Geometry parameters
  G4double worldSize  = 1000*1*angstrom;

  if ( !fpDefaultMaterial )
  {
    G4cerr << "Cannot retrieve materials already defined. " << G4endl;
    G4cerr << "Exiting application " << G4endl;
    exit(1);
  }

  //
  // World
  //
  G4VSolid* worldS
  = new G4Box("World",           // its name
              worldSize/2, worldSize/2, worldSize/2); // its size

  G4LogicalVolume*
  worldLV
  = new G4LogicalVolume(
      worldS,           // its solid
      fpDefaultMaterial,  // its material
      "World");         // its name

  G4VisAttributes *MyVisAtt_ZZ = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));
  MyVisAtt_ZZ ->SetVisibility (false);
  worldLV->SetVisAttributes(MyVisAtt_ZZ);

  G4VPhysicalVolume*
  worldPV
  = new G4PVPlacement(
      0,                // no rotation
      G4ThreeVector(),  // at (0,0,0)
      worldLV,          // its logical volume
      "World",          // its name
      0,                // its mother  volume
      false,            // no boolean operation
      0,                // copy number
      true);            // checking overlaps forced to YES

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes(G4String filename,
    unsigned short int option)
{
  //Clean old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Define Materials
  //
  ConstructMaterials();

  //Reconstruct world not to superimpose on geometries
  G4VPhysicalVolume* worldPV;
  G4LogicalVolume* worldLV;
  worldPV=ConstructWorld();
  worldLV=worldPV->GetLogicalVolume();

  //List of molecules inside pdb file separated by TER keyword
  fpMoleculeList=NULL;

  //'fPDBlib.load' is currently done for each 'DefineVolumes' call.
  int verbosity=0;  //To be implemented
  unsigned short int isDNA;
  fpMoleculeList = fPDBlib.Load(filename,isDNA,verbosity);



  if (fpMoleculeList!=NULL)
    {
    fPdbFileStatus=1;
    }

  if (option==1)
  {
    AtomisticView( worldLV,fpMoleculeList,1.0);
  }

      else
        if (option==10)
        {
          DrawBoundingVolume( worldLV,fpMoleculeList);
        }
        else
          if (option==11)
          {
            AtomisticView( worldLV,fpMoleculeList,1.0);
            DrawBoundingVolume( worldLV,fpMoleculeList);
          }

  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PDBlib DetectorConstruction::GetPDBlib()
{
  return fPDBlib;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule *DetectorConstruction::GetMoleculeList()
{
  return fpMoleculeList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////
///////////////// BEGIN atomistic representation
//
void DetectorConstruction::AtomisticView(G4LogicalVolume* worldLV,
    Molecule *moleculeListTemp, double atomSizeFactor)
{
  CheckMaterials();

  // All solids are defined for different color attributes :
  G4double sphereSize  = atomSizeFactor*1*angstrom;
  G4VSolid* atomS_H = new G4Orb("Sphere", sphereSize*1.2);
  G4VSolid* atomS_C = new G4Orb("Sphere", sphereSize*1.7);
  G4VSolid* atomS_O = new G4Orb("Sphere", sphereSize*1.52);
  G4VSolid* atomS_N = new G4Orb("Sphere", sphereSize*1.55);
  G4VSolid* atomS_S = new G4Orb("Sphere", sphereSize*1.8);
  G4VSolid* atomS_P = new G4Orb("Sphere", sphereSize*1.8);
  G4VSolid* atomS_X = new G4Orb("Sphere", sphereSize);  //Undefined

  //Logical volumes and color table for atoms :
  G4LogicalVolume* atomLV_H = new G4LogicalVolume(
      atomS_H, fpWaterMaterial,  "atomLV_H");  // its solid, material, name
  G4VisAttributes * MyVisAtt_H = new G4VisAttributes(
      G4Colour(G4Colour::White()));
  MyVisAtt_H->SetForceSolid(true);
  atomLV_H->SetVisAttributes(MyVisAtt_H);

  G4LogicalVolume* atomLV_C = new G4LogicalVolume(
      atomS_C, fpWaterMaterial, "atomLV_C");  // its solid, material, name
  G4VisAttributes * MyVisAtt_C = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));//Black in CPK, but bad rendering
  MyVisAtt_C->SetForceSolid(true);
  atomLV_C->SetVisAttributes(MyVisAtt_C);

  G4LogicalVolume* atomLV_O = new G4LogicalVolume(
      atomS_O, fpWaterMaterial, "atomLV_O"); // its solid, material, name
  G4VisAttributes * MyVisAtt_O = new G4VisAttributes(
      G4Colour(G4Colour::Red()));
  MyVisAtt_O->SetForceSolid(true);
  atomLV_O->SetVisAttributes(MyVisAtt_O);

  G4LogicalVolume* atomLV_N = new G4LogicalVolume(
      atomS_N, fpWaterMaterial, "atomLV_N"); // its solid, material, name
  G4VisAttributes * MyVisAtt_N = new G4VisAttributes(
      G4Colour(G4Colour(0.,0.,0.5)));//Dark blue
  MyVisAtt_N->SetForceSolid(true);
  atomLV_N->SetVisAttributes(MyVisAtt_N);

  G4LogicalVolume* atomLV_S = new G4LogicalVolume(
      atomS_S, fpWaterMaterial, "atomLV_S"); // its solid, material, name
  G4VisAttributes * MyVisAtt_S = new G4VisAttributes(G4Colour(
      G4Colour::Yellow()));
  MyVisAtt_S->SetForceSolid(true);
  atomLV_S->SetVisAttributes(MyVisAtt_S);

  G4LogicalVolume* atomLV_P = new G4LogicalVolume(
      atomS_P, fpWaterMaterial, "atomLV_P"); // its solid, material, name
  G4VisAttributes * MyVisAtt_P = new G4VisAttributes(
      G4Colour(G4Colour(1.0,0.5,0.)));//Orange
  MyVisAtt_P->SetForceSolid(true);
  atomLV_P->SetVisAttributes(MyVisAtt_P);

  G4LogicalVolume* atomLV_X = new G4LogicalVolume(
      atomS_X, fpWaterMaterial, "atomLV_X"); // its solid, material, name
  G4VisAttributes * MyVisAtt_X = new G4VisAttributes(
      G4Colour(G4Colour(1.0,0.75,0.8)));//Pink, other elements in CKP
  MyVisAtt_X->SetForceSolid(true);
  atomLV_X->SetVisAttributes(MyVisAtt_X);

  //Placement and physical volume construction from memory
  Residue *residueListTemp;
  Atom *AtomTemp;

  int nbAtomTot=0;  //Number total of atoms
  int nbAtomH=0, nbAtomC=0, nbAtomO=0, nbAtomN=0, nbAtomS=0, nbAtomP=0;
  int nbAtomX=0;
  int k=0;

  while (moleculeListTemp)
  {
    residueListTemp = moleculeListTemp->GetFirst();

    k++;
    int j=0;
    while (residueListTemp)
    {
      AtomTemp=residueListTemp->GetFirst();
      j++;

      int startFrom=0;
      int upTo=residueListTemp->fNbAtom; //case Base+sugar+phosphat (all atoms)
      for (int i=0 ; i < (upTo + startFrom) ; i++) //Phosphat or Sugar or Base
      {

        if (AtomTemp->fElement.compare("H") == 0)
        {
          nbAtomH++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_H,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("C") == 0)
        {
          nbAtomC++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_C,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("O") == 0)
        {
          nbAtomO++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_O,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("N") == 0)
        {
          nbAtomN++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_N,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("S") == 0)
        {
          nbAtomS++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_S,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else if (AtomTemp->fElement.compare("P") == 0)
        {
          nbAtomP++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_P,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }
        else
        {
          nbAtomX++;
          new G4PVPlacement(0,
                            G4ThreeVector(AtomTemp->fX*1*angstrom,
                                          AtomTemp->fY*1*angstrom,
                                          AtomTemp->fZ*1*angstrom),
                                          atomLV_X,
                                          "atomP",
                                          worldLV,
                                          false,
                                          0,
                                          fCheckOverlaps);
        }

        nbAtomTot++;
        //}//End if (i>=startFrom)
        AtomTemp=AtomTemp->GetNext();
      }//end of for (  i=0 ; i < residueListTemp->nbAtom ; i++)

      residueListTemp=residueListTemp->GetNext();
    }//end of while (residueListTemp)

    moleculeListTemp=moleculeListTemp->GetNext();
  }//end of while (moleculeListTemp)

  G4cout << "**************** atomisticView(...) ****************" <<G4endl;
  G4cout << "Number of loaded chains = " << k <<G4endl;
  G4cout << "Number of Atoms      = " << nbAtomTot <<G4endl;
  G4cout << "Number of Hydrogens  = " << nbAtomH <<G4endl;
  G4cout << "Number of Carbons    = " << nbAtomC <<G4endl;
  G4cout << "Number of Oxygens    = " << nbAtomO <<G4endl;
  G4cout << "Number of Nitrogens  = " << nbAtomN <<G4endl;
  G4cout << "Number of Sulfurs    = " << nbAtomS <<G4endl;
  G4cout << "Number of Phosphorus = " << nbAtomP <<G4endl;
  G4cout << "Number of undifined atoms =" << nbAtomX <<G4endl<<G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////
///////////////// BEGIN draw a bounding volume
//
void DetectorConstruction::DrawBoundingVolume(G4LogicalVolume* worldLV,
    Molecule *moleculeListTemp)
{
  CheckMaterials();

  double dX,dY,dZ;//Dimensions for bounding volume
  double tX,tY,tZ;//Translation for bounding volume
  fPDBlib.ComputeBoundingVolumeParams(moleculeListTemp,
                                      dX, dY, dZ,
                                      tX, tY, tZ);

  //[Geometry] Build a box :
  G4VSolid* boundingS
  = new G4Box("Bounding", dX*1*angstrom, dY*1*angstrom, dZ*1*angstrom);

  G4LogicalVolume* boundingLV
  = new G4LogicalVolume(boundingS, fpWaterMaterial, "BoundingLV");

  G4RotationMatrix *pRot = new G4RotationMatrix();

  G4VisAttributes *MyVisAtt_ZZ = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));
  boundingLV->SetVisAttributes(MyVisAtt_ZZ);

  new G4PVPlacement(pRot,
                    G4ThreeVector(
                        tX*1*angstrom,
                        tY*1*angstrom,
                        tZ*1*angstrom),  // at (x,y,z)
                        boundingLV,
                        "boundingPV",
                        worldLV
                        ,false,
                        0,
                        fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::LoadPDBfile(G4String fileName)
{
  G4cout<<"Load PDB file : "<<fileName<<"."<<G4endl<<G4endl;
  fPdbFileName=fileName;
#ifdef G4MULTITHREADED
  G4MTRunManager::GetRunManager()->DefineWorldVolume(
      DefineVolumes(fPdbFileName,fChosenOption)
  );
  G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
  G4RunManager::GetRunManager()->DefineWorldVolume(
      DefineVolumes(fPdbFileName,fChosenOption)
  );
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildBoundingVolume()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
    G4cout<<"Build only world volume and bounding volume"
        " for computation."<<G4endl<<G4endl;
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,10)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,10)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawAtoms_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,1)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,1)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DrawAtomsWithBounding_()
{
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,11)
    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(fPdbFileName,11)
    );
#endif
  }
  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}



