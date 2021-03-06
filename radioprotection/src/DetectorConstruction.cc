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

#include "DetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"
#include "SensitiveDetector.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4AssemblyVolume.hh"

using namespace std;

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager):
	solidWP(0),
	solidTarget(0),
	logicTreatmentRoom(0),
	logicWP(0), 
	logicTarget(0), 
	logicScintStart_1(0), 
	logicScintVeto_1(0), 
	logicScintVeto_2(0),
	logicLysoAlShell(0), 
	logicLyso(0), 
	logicNaIAlShell(0), 
	logicNaI(0), 
	physicalTreatmentRoom(0), 
	physicalWP(0), 
	physicalTarget(0),
	physicalScintStart_1(0), 
	physicalScintVeto_1(0), 
	physicalScintVeto_2(0),
	physicalLysoAlShell(0), 
	physicalLyso(0), 
	physicalNaIAlShell(0), 
	physicalNaI(0) 
{
	analysis = analysis_manager;
	DetectorMessenger = new DetectorConstructionMessenger(this);
}

DetectorConstruction::~DetectorConstruction(){

}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

	SetDefaultDimensions();

	// -----------------------------
	// Treatment room - World volume
	//------------------------------
	// Treatment room sizes
	const G4double worldX = 100.0 *cm;
	const G4double worldY = 100.0 *cm;
	const G4double worldZ = 300.0 *cm;
	G4bool isotopes = false;

	G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);

	G4Box* treatmentRoom = new G4Box("TreatmentRoom",
			worldX*0.5,
			worldY*0.5,
			worldZ*0.5);

	logicTreatmentRoom = new G4LogicalVolume(treatmentRoom,
			airNist,
			"logicTreatmentRoom",
			0,0,0);

	physicalTreatmentRoom = new G4PVPlacement(0,
			G4ThreeVector(),
			"physicalTreatmentRoom",
			logicTreatmentRoom,
			0,
			false,
			0);


	// The treatment room is invisible in the Visualisation
	logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::GetInvisible());


	NaI();
	LaBr(); //substite of LYSO                
	//LYSO();
	Target();
	WaterPhantom();
	ScintStart_1();
	ScintVeto_1();
	ScintVeto_2();


	return physicalTreatmentRoom; 

}

void DetectorConstruction::WaterPhantom()
{
	/***** WATER PHANTOM *****/
	G4bool isotopes = false;
	//G4Material* WaterNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);

	solidWP = new G4Box("solidWP",
			WP_HLX,
			WP_HLY,
			WP_HLZ);

	logicWP = new G4LogicalVolume(solidWP,
			WP_mat,
			"logicWP",
			0,0,0);

	physicalWP = new G4PVPlacement(0,                     //no rotation
			WP_Trans,       //translation
			logicWP,            //its logical volume
			"physicalWP",               //its name
			logicTreatmentRoom,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking
}

void DetectorConstruction::Target()
{
	solidTarget = new G4Tubs("solidTarget",                       //its name
			Target_RMin, 
			Target_RMax, 
			Target_HLZ,
			0,
			360*deg);     //its size

	logicTarget = new G4LogicalVolume(solidTarget,          //its solid
			Target_mat,           //its material
			"logicTarget");            //its name

	physicalTarget =  new G4PVPlacement(G4Transform3D(TargetRot,  Target_Trans),       //la faccia in z del WP e' nell'origine
			logicTarget,            //its logical volume
			"physicalTarget",               //its name
			logicTreatmentRoom,   //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking

}

/////////////////////////////////////////////////////////////////////////////
void DetectorConstruction::LYSO()
{
	// ------------------//
	// WATER PHANTOM   //
	//-------------------//

	G4Box* solidLysoAlShell = new G4Box("solidLysoAlShell",                       //its name
			LysoAlShell_HLX, 
			LysoAlShell_HLY, 
			LysoAlShell_HLZ);     //its size

	logicLysoAlShell = new G4LogicalVolume(solidLysoAlShell,          //its solid
			LysoAlShell_mat,           //its material
			"LysoAlShellLV");            //its name

	physicalLysoAlShell =  new G4PVPlacement( G4Transform3D(YRot, Lyso_Trans),      //la faccia in z del WP e' nell'origine
			logicLysoAlShell,            //its logical volume
			"physicalLysoAlShell",               //its name
			logicTreatmentRoom,   //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking
	//logicLysoAlShell -> SetVisAttributes(skyBlue);


	G4Box* solidLyso = new G4Box("solidLyso",                       //its name
			Lyso_HLX, 
			Lyso_HLY, 
			Lyso_HLZ);     //its size

	logicLyso = new G4LogicalVolume(solidLyso,          //its solid
			Lyso_mat,           //its material
			"logicLyso");            //its name

	physicalLyso =  new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(0,0,0),       //la faccia in z del WP e' nell'origine
			logicLyso,            //its logical volume
			"physicalLyso",               //its name
			logicLysoAlShell,   //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking
	//logicLyso -> SetVisAttributes(red);


}


//////////////////////////////////////////////////////////////
void DetectorConstruction::NaI()
{
	// ------------------//
	// WATER PHANTOM	 //
	//-------------------//

	//Lyso_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

	G4Tubs* solidNaIAlShell = new G4Tubs("solidNaIAlShell",                       //its name
			NaIAlShell_RMin, 
			NaIAlShell_RMax, 
			NaIAlShell_HLZ,
			0,
			360*deg);     //its size

	logicNaIAlShell = new G4LogicalVolume(solidNaIAlShell,          //its solid
			NaIAlShell_mat,           //its material
			"logicNaIAlShell");            //its name

	physicalNaIAlShell =  new G4PVPlacement( G4Transform3D(YRot, NaI_Trans),       //la faccia in z del WP e' nell'origine
			logicNaIAlShell,            //its logical volume
			"physicalNaIAlShell",               //its name
			logicTreatmentRoom,   //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking


	G4Tubs* solidNaI = new G4Tubs("solidNaI",                       //its name
			NaI_RMin, 
			NaI_RMax, 
			NaI_HLZ,
			0,
			360*deg);     //its size

	logicNaI = new G4LogicalVolume(solidNaI,          //its solid
			NaI_mat,           //its material
			"logicNaI");            //its name

	physicalNaI =  new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(0,0,0),       //la faccia in z del WP e' nell'origine
			logicNaI,            //its logical volume
			"physicalNaI",               //its name
			logicNaIAlShell,   //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking


}

//////////////////////////////////////////////////////////////
void DetectorConstruction::LaBr()
{
	// ------------------//
	// WATER PHANTOM	 //
	//-------------------//

	//Lyso_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

	G4Tubs* solidLaBrAlShell = new G4Tubs("solidLaBrAlShell",                       //its name
			LaBrAlShell_RMin, 
			LaBrAlShell_RMax, 
			LaBrAlShell_HLZ,
			0,
			360*deg);     //its size

	logicLaBrAlShell = new G4LogicalVolume(solidLaBrAlShell,          //its solid
			LaBrAlShell_mat,           //its material
			"logicLaBrAlShell");            //its name

	physicalLaBrAlShell =  new G4PVPlacement( G4Transform3D(YRot, LaBr_Trans),       //la faccia in z del WP e' nell'origine
			logicLaBrAlShell,            //its logical volume
			"physicalLaBrAlShell",               //its name
			logicTreatmentRoom,   //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking


	G4Tubs* solidLaBr = new G4Tubs("solidLaBr",                       //its name
			LaBr_RMin, 
			LaBr_RMax, 
			LaBr_HLZ,
			0,
			360*deg);     //its size

	logicLaBr = new G4LogicalVolume(solidLaBr,          //its solid
			LaBr_mat,           //its material
			"logicLaBr");            //its name

	physicalLaBr =  new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(0,0,0),       //la faccia in z del WP e' nell'origine
			logicLaBr,            //its logical volume
			"physicalLaBr",               //its name
			logicLaBrAlShell,   //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking


}


////////////////////////////////////////////////////////////////////////
void DetectorConstruction::ScintStart_1()
{

	G4Box* solidScintStart_1 = new G4Box("solidScintStart_1",
			ScintStart_1_HLX,
			ScintStart_1_HLY,
			ScintStart_1_HLZ);

	logicScintStart_1 = new G4LogicalVolume(solidScintStart_1,
			BC400,
			"logicScintStart_1",
			0,0,0);

	physicalScintStart_1 = new G4PVPlacement(0,                     //no rotation
			ScintStart_1_Trans,       //translation
			logicScintStart_1,            //its logical volume
			"physicalScintStart_1",               //its name
			logicTreatmentRoom,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking


}

////////////////////////////////////////////////////////////////////////
void DetectorConstruction::ScintVeto_1()
{

	G4Box* solidScintVeto_1 = new G4Box("solidScintVeto_1",
			ScintVeto_1_HLX,
			ScintVeto_1_HLY,
			ScintVeto_1_HLZ);

	logicScintVeto_1 = new G4LogicalVolume(solidScintVeto_1,
			BC400,
			"logicScintVeto_1",
			0,0,0);

	physicalScintVeto_1 = new G4PVPlacement(G4Transform3D(YRot, ScintVeto_1_Trans),      //translation
			logicScintVeto_1,            //its logical volume
			"physicalScintVeto_1",               //its name
			logicTreatmentRoom,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking

}

////////////////////////////////////////////////////////////////////////
void DetectorConstruction::ScintVeto_2()
{

	G4Box* solidScintVeto_2 = new G4Box("solidScintVeto_2",
			ScintVeto_2_HLX,
			ScintVeto_2_HLY,
			ScintVeto_2_HLZ);

	logicScintVeto_2 = new G4LogicalVolume(solidScintVeto_2,
			BC400,
			"logicScintVeto_2",
			0,0,0);

	physicalScintVeto_2 = new G4PVPlacement(G4Transform3D(YRot, ScintVeto_2_Trans),      //translation
			logicScintVeto_2,            //its logical volume
			"physicalScintVeto_2",               //its name
			logicTreatmentRoom,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking

}
////////////////////////////////////////////////////////////////////////
void DetectorConstruction::SetDefaultDimensions()
{

	// ELEMENTS
	G4bool isotopes = false;
	G4Material* Al_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

	ScintStart_1_HLX = 2.5*cm;
	ScintStart_1_HLY = 2.5*cm;
	ScintStart_1_HLZ = 1.5*mm;
	ScintStart_1_Trans = G4ThreeVector(0,0,77*cm);


	WP_HLX = 5.*cm;
	WP_HLY = 5.*cm;
	WP_HLZ = 1*cm;
	WP_Trans = G4ThreeVector(0,0,40*cm);
	WP_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

	Target_RMax = 0.5*cm; //Y89Default
	Target_RMin = 0.*mm;
	Target_HLZ = 5*mm;  //1cm alt0
	//TargetRot.rotateX(90*deg);
	TargetRot.rotateX(0*deg);
	Target_Trans = G4ThreeVector(0,0,120*cm); //the NaI entrance is at 89 cm from WP lateral edge
	Target_mat = MaterialWithSingleIsotope("YttriumRod","YRod", 4.47*g/cm3, 39,  89);

	G4Element* el_Cu =  new G4Element("El_Copper","Cu", 29,  63.54*g/mole);
	Cu_mat = new G4Material("Copper_mat", 8.92*g/cm3, 1);
	Cu_mat -> AddElement(el_Cu, 1);


	ScintVeto_1_HLX = 2.5*cm;
	ScintVeto_1_HLY = 2.5*cm;
	ScintVeto_1_HLZ = 2.5*mm;
	ScintVeto_1_Trans = Target_Trans + G4ThreeVector(15.*cm,0,0);

	ScintVeto_2_HLX = 2.5*cm;
	ScintVeto_2_HLY = 2.5*cm;
	ScintVeto_2_HLZ = 2.5*mm;
	ScintVeto_2_Trans = Target_Trans - G4ThreeVector(15.*cm,0,0);


	//LysoCrystal
	LysoAlShell_HLX = 22.*mm;
	LysoAlShell_HLY = 22.*mm;
	LysoAlShell_HLZ = 41.*mm;
	Lyso_HLX = 20.*mm;
	Lyso_HLY = 20.*mm;
	Lyso_HLZ = 40.*mm;
	YRot.rotateY(90*deg);
	Lyso_Trans = ScintVeto_1_Trans + G4ThreeVector(3.5*cm+LysoAlShell_HLZ,0,0); //the lyso entrance is at 89 cm from WP lateral edge

	LysoAlShell_mat = Al_mat;
	G4Element* el_Lu =  new G4Element("Lutezio","Lu", 71,  174.967*g/mole);
	G4Element* el_Y =  new G4Element("Yttrium","Y", 39,  88.90585*g/mole);
	G4Element* el_Si =  new G4Element("Silicon","Si", 14,  28.0855*g/mole);
	G4Element* el_O =  new G4Element("Oxygen","O", 8,  15.999*g/mole);
	Lyso_mat = new G4Material("LYSO", 7.2*g/cm3, 4);
	Lyso_mat -> AddElement(el_Lu, 0.680);
	Lyso_mat -> AddElement(el_Y, 0.044);
	Lyso_mat -> AddElement(el_Si, 0.072);
	Lyso_mat -> AddElement(el_O, 0.204);

	//NaI
	NaIAlShell_RMax = 25.9*mm;
	NaIAlShell_RMin = 0.*mm;
	NaIAlShell_HLZ = 27*mm;
	NaI_RMax = 25.5*mm;
	NaI_RMin = 0.*mm;
	NaI_HLZ = 25.5*mm;
	NaI_Trans = ScintVeto_2_Trans - G4ThreeVector(3.5*cm+NaIAlShell_HLZ,0,0); //the NaI entrance is at 89 cm from WP lateral edge

	NaIAlShell_mat = Al_mat;
	NaI_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");

	//LaBr
	G4Element* elLa = new G4Element("Lanthanum", "La", 57, 138.905*g/mole);
        G4Element* elBr = new G4Element("Bromine","Br", 35, 79.904*g/mole);
	LaBr_mat = new G4Material("LaBr3",5.08*g/cm3, 2);
        LaBr_mat->AddElement(elLa, 1);
        LaBr_mat->AddElement(elBr, 3);
	
	LaBrAlShell_RMax = 43.2*mm; //3mm
	LaBrAlShell_RMin = 0.*mm;
	LaBrAlShell_HLZ = 44.9*mm;
	LaBr_RMax = 40.2*mm;
	LaBr_RMin = 0.*mm;
	LaBr_HLZ = 41.9*mm;
	LaBr_Trans = ScintVeto_1_Trans + G4ThreeVector(3.5*cm+LaBrAlShell_HLZ,0,0); //the LaBr entrance is at 89 cm from WP lateral edge

	LaBrAlShell_mat = Al_mat;
	LaBr_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");

	// COMPOUND
	G4Element* elH = new G4Element("Hydrogen", "H", 1., 1.01*g/mole);
	G4Element* elC = new G4Element("Oxygen","C", 6., 12.0107*g/mole);

	BC400 = new G4Material("BC400",1.032*g/cm3, 2);
	BC400->AddElement(elH, 0.08474);
	BC400->AddElement(elC, 0.91526);


	checkOverlaps = true;  //<---------------------
}



////////////////////////////////////////////////////////////////////////
G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
		G4String symbol, G4double density, G4int Z, G4int A)
{
	// define a material from an isotope
	//
	G4int ncomponents;
	G4double abundance, massfraction;

	G4Isotope* isotope = new G4Isotope(symbol, Z, A);

	G4Element* element  = new G4Element(name, symbol, ncomponents=1);
	element->AddIsotope(isotope, abundance= 100.*perCent);

	G4Material* material = new G4Material(name, density, ncomponents=1);
	material->AddElement(element, massfraction=100.*perCent);

	return material;
}


void DetectorConstruction::SetAbsorberSize(G4ThreeVector size)
{
	solidWP ->SetXHalfLength(size.getX());
	solidWP ->SetYHalfLength(size.getY());
	solidWP ->SetZHalfLength(size.getZ());
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
	G4cout<<"The x y z size of the WP  (mm):\t"
		<< ((solidWP -> GetXHalfLength())*2) << " " << ((solidWP -> GetYHalfLength())*2) << " " << ((solidWP -> GetZHalfLength())*2)<< G4endl;

}
void DetectorConstruction::SetAbsorberMaterial(G4String material)
{
	if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
	{
		if (pttoMaterial)
		{
			WP_mat  = pttoMaterial;
			logicWP -> SetMaterial(WP_mat);
			G4cout << "The material of the Absorber has been changed to " << material << G4endl;
			G4RunManager::GetRunManager()->GeometryHasBeenModified();
		}
	}
	else
	{
		G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
			" table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
	}
}
////////////////////////////////////////////////////////////
void DetectorConstruction::SetTargetMaterial(G4String material)
{
	if(material == "Copper")
	{
		Target_mat = Cu_mat;
		logicTarget ->SetMaterial(Target_mat);
		G4cout << "The material of the Absorber has been changed to " << material << G4endl;
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
	}
	else if(material == "Yttrium")
	{
		logicTarget ->SetMaterial(Target_mat);
		G4cout << "The material of the Absorber has been changed to " << material << G4endl;
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
	}
	else if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
	{
		if (pttoMaterial)
		{
			Target_mat  = pttoMaterial;
			logicTarget -> SetMaterial(Target_mat);
			G4cout << "The material of the Absorber has been changed to " << material << G4endl;
			G4RunManager::GetRunManager()->GeometryHasBeenModified();
		}
	}
	else
	{
		G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
			" table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
	}

}

void DetectorConstruction::SetTargetLength(G4double value)
{
	//solidTarget -> SetOuterRadius(value);
	solidTarget -> SetZHalfLength(value);
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

////////////////////////////////////////////////////////////////////////
void DetectorConstruction::ConstructSDandField()
{

	/*
	   LysoSD = new SensitiveDetector("LysoSD", "LysoHitsCollection", analysis);
	   G4SDManager::GetSDMpointer()->AddNewDetector(LysoSD); 
	   SetSensitiveDetector("logicLyso", LysoSD);

	   NaISD = new SensitiveDetector("NaIDetectorSD", "NaIHitsCollection", analysis);
	   G4SDManager::GetSDMpointer()->AddNewDetector(NaISD); 
	   SetSensitiveDetector("logicNaI", NaISD);

	   ScintStart_1SD = new SensitiveDetector("ScintStart_1SD", "ScintStart_1HitsCollection", analysis);
	   G4SDManager::GetSDMpointer()->AddNewDetector(ScintStart_1SD); 
	   SetSensitiveDetector("logicScintStart_1", ScintStart_1SD);

	   ScintVeto_1SD = new SensitiveDetector("ScintVeto_1SD", "ScintVeto_1HitsCollection", analysis);
	   G4SDManager::GetSDMpointer()->AddNewDetector(ScintVeto_1SD); 
	   SetSensitiveDetector("logicScintVeto_1", ScintVeto_1SD);

	   ScintVeto_2SD = new SensitiveDetector("ScintVeto_2SD", "ScintVeto_2HitsCollection", analysis);
	   G4SDManager::GetSDMpointer()->AddNewDetector(ScintVeto_2SD); 
	   SetSensitiveDetector("logicScintVeto_2", ScintVeto_2SD);

	 */


}
