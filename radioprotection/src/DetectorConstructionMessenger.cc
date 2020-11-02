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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy


#include "DetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

   DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* detector)
:TreatmentRoomGeom(detector)

{

    GeomDir = new G4UIdirectory("/changeGeometry/");
    GeomDir -> SetGuidance("Command to Change geometry");


    Target_Radius = new G4UIcmdWithADoubleAndUnit("/changeGeometry/TargetRadius", this);
    Target_Radius -> SetGuidance("Change the target radius");
    Target_Radius -> SetParameterName("TargetRadius", false);
    Target_Radius -> SetDefaultUnit("cm");
    Target_Radius -> SetUnitCandidates("nm um mm cm");
    Target_Radius -> AvailableForStates(G4State_Idle);

    Target_material = new G4UIcmdWithAString("/changeGeometry/TargetMaterial", this);
    Target_material -> SetGuidance("Change the target material");
    Target_material -> SetParameterName("TargetMaterial", false);
    Target_material -> SetDefaultValue("Yttrium");
    Target_material -> AvailableForStates(G4State_Idle);


    Absorber_material = new G4UIcmdWithAString("/changeGeometry/AbsorberMaterial", this);
    Absorber_material -> SetGuidance("Change the Absorber material");
    Absorber_material -> SetParameterName("AbsorberMaterial", false);
    Absorber_material -> SetDefaultValue("G4_WATER");
    Absorber_material -> AvailableForStates(G4State_Idle);

    // Change Phantom position
    Absorber_size = new G4UIcmdWith3VectorAndUnit("/changeGeometry/AbsorberSize", this);
    Absorber_size -> SetGuidance("Insert sizes X Y and Z\n   0 or negative values mean <<Don't change it!>>");
    Absorber_size -> SetParameterName("HLX", "HLY", "HLZ", false);
    Absorber_size -> SetDefaultUnit("mm");
    Absorber_size -> SetUnitCandidates("nm um mm cm");
    Absorber_size -> AvailableForStates(G4State_Idle);



}

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{ 
    delete Target_material;
    delete Target_Radius;
    delete Absorber_material;
    delete Absorber_size;
    delete GeomDir; 

}




void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
    if( command == Absorber_size )
    {
	    TreatmentRoomGeom -> SetAbsorberSize(Absorber_size -> GetNew3VectorValue(newValue));
    }
    else if( command == Absorber_material )
    {
                TreatmentRoomGeom -> SetAbsorberMaterial(newValue);
    }
    else if( command == Target_material)
    {
        TreatmentRoomGeom -> SetTargetMaterial(newValue);
    }
    else if( command == Target_Radius)
    {
        TreatmentRoomGeom -> SetTargetRadius(Target_Radius->GetNewDoubleValue(newValue));
    }


}

