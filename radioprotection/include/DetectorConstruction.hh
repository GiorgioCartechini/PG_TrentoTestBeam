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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#ifndef DetectorConstruction_H 
#define DetectorConstruction_H 1
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Trd.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"


class G4VPhysicalVolume;
class DetectorMessenger;
class G4LogicalVolume;
class SensitiveDetector;
class DetectorConstructionMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		DetectorConstruction(AnalysisManager* analysis);
		~DetectorConstruction();

		G4VPhysicalVolume* Construct();

		void WaterPhantom();
		void Target();
		void LYSO();
		void NaI();
		void LaBr();
		void ScintStart_1();
		void ScintVeto_1();
		void ScintVeto_2();

		G4Material* MaterialWithSingleIsotope( G4String name, G4String symbol, G4double density, G4int Z, G4int A);

		void ConstructSDandField();

		void SetAbsorberSize(G4ThreeVector size);
		void SetAbsorberMaterial(G4String material);
		void SetTargetMaterial(G4String material);
		void SetTargetLength(G4double radius);

	private:
		AnalysisManager* analysis;

		DetectorConstructionMessenger* DetectorMessenger;

		void SetDefaultDimensions();


		G4Box* solidWP;
		G4Tubs* solidTarget;

		G4LogicalVolume* logicTreatmentRoom;
		G4LogicalVolume* logicWP;
		G4LogicalVolume* logicTarget;
		G4LogicalVolume* logicScintStart_1;
		G4LogicalVolume* logicScintVeto_1;
		G4LogicalVolume* logicScintVeto_2;
		G4LogicalVolume* logicLysoAlShell;
		G4LogicalVolume* logicLyso;
		G4LogicalVolume* logicNaIAlShell;
		G4LogicalVolume* logicNaI;
		G4LogicalVolume* logicLaBrAlShell;
		G4LogicalVolume* logicLaBr;

		G4VPhysicalVolume* physicalTreatmentRoom;
		G4VPhysicalVolume* physicalWP;
		G4VPhysicalVolume* physicalTarget;
		G4VPhysicalVolume* physicalScintStart_1;
		G4VPhysicalVolume* physicalScintVeto_1;
		G4VPhysicalVolume* physicalScintVeto_2;
		G4VPhysicalVolume* physicalLysoAlShell;
		G4VPhysicalVolume* physicalLyso;
		G4VPhysicalVolume* physicalNaIAlShell;
		G4VPhysicalVolume* physicalNaI;
		G4VPhysicalVolume* physicalLaBrAlShell;
		G4VPhysicalVolume* physicalLaBr;




		G4double ScintStart_1_HLX;
		G4double ScintStart_1_HLY;
		G4double ScintStart_1_HLZ;

		G4double ScintVeto_1_HLX;
		G4double ScintVeto_1_HLY;
		G4double ScintVeto_1_HLZ;

		G4double ScintVeto_2_HLX;
		G4double ScintVeto_2_HLY;
		G4double ScintVeto_2_HLZ;

		G4double WP_HLX;
		G4double WP_HLY;
		G4double WP_HLZ;

		G4double Target_RMax;
		G4double Target_RMin;
		G4double Target_HLZ;

		G4double LysoAlShell_HLX;
		G4double LysoAlShell_HLY;
		G4double LysoAlShell_HLZ;

		G4double Lyso_HLX;
		G4double Lyso_HLY;
		G4double Lyso_HLZ;

		G4double NaIAlShell_RMax;
		G4double NaIAlShell_RMin;
		G4double NaIAlShell_HLZ;
		G4double NaI_RMax;
		G4double NaI_RMin;
		G4double NaI_HLZ;

		G4double LaBrAlShell_RMax;
		G4double LaBrAlShell_RMin;
		G4double LaBrAlShell_HLZ;
		G4double LaBr_RMax;
		G4double LaBr_RMin;
		G4double LaBr_HLZ;


		//baf MATERIALS
		G4Material* BC400;
		G4Material* Target_mat;
		G4Material* WP_mat;
		G4Material* Lyso_mat;
		G4Material* LysoAlShell_mat;
		G4Material* NaIAlShell_mat;
		G4Material* LaBrAlShell_mat;
		G4Material* NaI_mat; 
		G4Material* LaBr_mat;
		G4Material* Cu_mat;


		G4RotationMatrix YRot;
		G4RotationMatrix TargetRot;
		G4ThreeVector ScintStart_1_Trans;
		G4ThreeVector ScintVeto_1_Trans;
		G4ThreeVector ScintVeto_2_Trans;
		G4ThreeVector WP_Trans;
		G4ThreeVector Lyso_Trans;
		G4ThreeVector NaI_Trans;
		G4ThreeVector LaBr_Trans;
		G4ThreeVector Target_Trans;

		SensitiveDetector* NaISD;
		SensitiveDetector* LysoSD;
		SensitiveDetector* ScintStart_1SD;
		SensitiveDetector* ScintVeto_1SD;
		SensitiveDetector* ScintVeto_2SD;

		G4bool checkOverlaps;
};
#endif
