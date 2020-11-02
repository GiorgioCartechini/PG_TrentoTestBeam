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
//  code based on the basic example B2
//
#include "SensitiveDetector.hh"
#include "AnalysisManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4StepPoint.hh"

#include "EventAction.hh"

SensitiveDetector::SensitiveDetector(const G4String& name,
                         const G4String& hitsCollectionName, AnalysisManager* analysis_manager) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL),
   fnewHit(0),
    fTime(-10),
    fEkin(-10),
    fTotalEdep(-10),
    fPDG(-10)
{
  collectionName.insert(hitsCollectionName);
  analysis = analysis_manager;
  
}

SensitiveDetector::~SensitiveDetector() 
{}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new SensitiveDetectorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce 

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  fnewHit = new SensitiveDetectorHit();
  //fHitsCollection->insert(new SensitiveDetectorHit());

  //fnewHit = (*fHitsCollection)[0];


  fTotalEdep = 0.;

}


G4bool SensitiveDetector::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory* /*history*/)
{  
  G4StepPoint* Prestep = aStep->GetPreStepPoint();
  //G4StepPoint* Poststep = aStep->GetPostStepPoint();

  G4String fvolumeName =  Prestep->GetTouchableHandle() ->GetVolume() -> GetName();

    
    if(Prestep ->GetStepStatus() == fGeomBoundary)
    {
      //G4cout << fvolumeName << G4endl;
      fEkin = Prestep->GetKineticEnergy();
      fTime = aStep  ->  GetTrack() -> GetGlobalTime();
      fPName =  aStep  ->  GetTrack() -> GetDefinition() -> GetParticleName();
      fPDG = aStep  ->  GetTrack() -> GetDynamicParticle() -> GetPDGcode();
      //G4cout <<volumeName <<"\t" << fPDG<<"\t" << fPName << "\t" << fTime/ns << G4endl;
    }

    fTotalEdep += aStep -> GetTotalEnergyDeposit();
	
	
	return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{

  fnewHit -> SetPDG(fPDG);
  fnewHit -> SetEkin(fEkin);
  fnewHit -> SetTime(fTime);
  fnewHit -> SetEdep(fTotalEdep);
  fHitsCollection->insert(fnewHit);

  /*if(SensitiveDetectorName == "LysoSD")
    analysis -> FillLyso(fnewHit -> GetPDG(), fnewHit -> GetEkin()/MeV, fnewHit -> GetEdep()/MeV, fnewHit -> GetTime()/ns);
  else if(SensitiveDetectorName == "NaISD")
    analysis -> FillNaI(fnewHit -> GetPDG(), fnewHit -> GetEkin()/MeV, fnewHit -> GetEdep()/MeV, fnewHit -> GetTime()/ns);
  else if(SensitiveDetectorName == "ScintStart_1SD")
    analysis -> FillScintStart1(fnewHit -> GetPDG(), fnewHit -> GetEkin()/MeV, fnewHit -> GetEdep()/MeV, fnewHit -> GetTime()/ns);
  else if(SensitiveDetectorName == "ScintVeto_1SD")
    analysis -> FillScintVeto1(fnewHit -> GetPDG(), fnewHit -> GetEkin()/MeV, fnewHit -> GetEdep()/MeV, fnewHit -> GetTime()/ns);
  else if(SensitiveDetectorName == "ScintVeto_2SD")
    analysis -> FillScintVeto2(fnewHit -> GetPDG(), fnewHit -> GetEkin()/MeV, fnewHit -> GetEdep()/MeV, fnewHit -> GetTime()/ns);*/

  //G4cout << SensitiveDetectorName <<" "<<fHitsCollection->entries() << G4endl;

}
