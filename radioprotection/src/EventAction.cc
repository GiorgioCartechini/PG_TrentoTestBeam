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
// 
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "SensitiveDetector.hh"
#include "SensitiveDetectorHit.hh"
#include "AnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(AnalysisManager* analysis)
	: G4UserEventAction(),
	fLysoHCID(-1),
	fNaIHCID(-1),
	fScintStart1HCID(-1),
	fScintVeto1HCID(-1),
	fScintVeto2HCID(-1)
{
	analysisMan = analysis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDetectorHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
		const G4Event* event) const
{
	auto hitsCollection 
		= static_cast<SensitiveDetectorHitsCollection*>(
				event->GetHCofThisEvent()->GetHC(hcID));

	if ( ! hitsCollection ) {
		G4ExceptionDescription msg;
		msg << "Cannot access hitsCollection ID " << hcID; 
		G4Exception("EventAction::GetHitsCollection()",
				"MyCode0003", FatalException, msg);
	}         

	return hitsCollection;
}    


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{

	fPDGScStart1 = -10;
	fVolScStart1 = -10;
	fEkinScStart1 =-10*MeV;
	fEdepScStart1 = 0.*MeV;
	fTimeScStart1 = -10*ns;

	fPDGScVeto1 = -10;
	fVolScVeto1 = -10;
	fEkinScVeto1 = -10*MeV;
	fEdepScVeto1 = 0.*MeV;
	fTimeScVeto1 = -10*ns;

	fPDGScVeto2 = -10;
	fVolScVeto2 = -10;
	fEkinScVeto2 = -10*MeV;
	fEdepScVeto2 = 0.*MeV;
	fTimeScVeto2 = -10*ns;

	fPDGNaI = -10;
	fVolNaI = -10;
	fEkinNaI =-10*MeV;
	fEdepNaI = 0.*MeV;
	fTimeNaI = -10*ns;

	fPDGLyso = -10;
	fVolLyso = -10;
	fEkinLyso =-10*MeV;
	fEdepLyso = 0.*MeV;
	fTimeLyso = -10*ns;

	fPDGLaBr = -10;
	fVolLaBr = -10;
	fEkinLaBr =-10*MeV;
	fEdepLaBr = 0.*MeV;
	fTimeLaBr = -10*ns;


	fPDGLyso = -10;
	fVolTgt = -10;
	fEkinLyso =-10*MeV;
	fEdepLyso = 0.*MeV;
	fTimeLyso = -10*ns;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{  



	//G4cout << GetPDGScStart1() <<"\t" << GetEkinScStart1()/MeV <<"\t" << GetEdepScStart1()/MeV << G4endl;
	//G4cout << GetPDGScVeto1() <<"\t" << GetEkinScVeto1()/MeV <<"\t" << GetEdepScVeto1()/MeV << G4endl;
	//G4cout << GetPDGScVeto2() <<"\t" << GetEkinScVeto2()/MeV <<"\t" << GetEdepScVeto2()/MeV << G4endl;
	//G4cout << GetPDGNaI() <<"\t" << GetEkinNaI()/MeV <<"\t" << GetEdepNaI()/MeV << G4endl;
	//G4cout << GetPDGLyso() <<"\t" << GetEkinLyso()/MeV <<"\t" << GetEdepLyso()/MeV << G4endl;

	//analysisMan -> FillScintStart1(GetPDGScStart1(), GetVertexVolScStart1(), GetEkinScStart1()/MeV);
	//analysisMan -> FillScintVeto1(GetPDGScVeto1(), GetVertexVolScVeto1(), GetEkinScVeto1()/MeV);
	//analysisMan -> FillScintVeto2(GetPDGScVeto2(), GetVertexVolScVeto2(), GetEkinScVeto2()/MeV);
	//analysisMan -> FillNaI(GetPDGNaI(), GetVertexVolNaI(), GetEkinNaI()/MeV);
	//analysisMan -> FillLyso(GetPDGLyso(), GetVertexVolLyso(), GetEkinLyso()/MeV);
	analysisMan -> FillNaIEdep(GetEdepNaI()/MeV);
	analysisMan -> FillLysoEdep(GetEdepLyso()/MeV);
	analysisMan -> FillLaBrEdep(GetEdepLaBr()/MeV);
	/*
	   if ( fLysoHCID == -1 )
	   fLysoHCID = G4SDManager::GetSDMpointer()->GetCollectionID("LysoHitsCollection");

	   if ( fNaIHCID == -1 )
	   fNaIHCID = G4SDManager::GetSDMpointer()->GetCollectionID("NaIHitsCollection");

	   if ( fScintStart1HCID == -1 )
	   fScintStart1HCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintStart_1HitsCollection");

	   if ( fScintVeto1HCID == -1 )
	   fScintVeto1HCID  = G4SDManager::GetSDMpointer()->GetCollectionID("ScintVeto_1HitsCollection");

	   if ( fScintVeto2HCID == -1 )
	   fScintVeto2HCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintVeto_2HitsCollection");

	// Get hits collections
	auto lysoHC = GetHitsCollection(fLysoHCID, event);
	auto NaIHC = GetHitsCollection(fNaIHCID, event);
	auto ScintStart1HC = GetHitsCollection(fScintStart1HCID, event);
	auto ScintVeto1HC = GetHitsCollection(fScintVeto1HCID, event);
	auto ScintVeto2HC = GetHitsCollection(fScintVeto2HCID, event);

	analysisMan -> FillLyso((*lysoHC)[0] -> GetPDG(), (*lysoHC)[0] -> GetEkin()/MeV, (*lysoHC)[0] -> GetEdep()/MeV, (*lysoHC)[0] -> GetTime()/ns);
	analysisMan -> FillNaI((*NaIHC)[0] -> GetPDG(), (*NaIHC)[0] -> GetEkin()/MeV, (*NaIHC)[0] -> GetEdep()/MeV, (*NaIHC)[0] -> GetTime()/ns);
	analysisMan -> FillScintStart1((*ScintStart1HC)[0] -> GetPDG(), (*ScintStart1HC)[0] -> GetEkin()/MeV, (*ScintStart1HC)[0] -> GetEdep()/MeV, (*ScintStart1HC)[0] -> GetTime()/ns);
	analysisMan -> FillScintVeto1((*ScintVeto1HC)[0] -> GetPDG(), (*ScintVeto1HC)[0] -> GetEkin()/MeV, (*ScintVeto1HC)[0] -> GetEdep()/MeV, (*ScintVeto1HC)[0] -> GetTime()/ns);
	analysisMan -> FillScintVeto2((*ScintVeto2HC)[0] -> GetPDG(), (*ScintVeto2HC)[0] -> GetEkin()/MeV, (*ScintVeto2HC)[0] -> GetEdep()/MeV, (*ScintVeto2HC)[0] -> GetTime()/ns);
	 */
	//G4cout <<(*ScintStart1HC)[0] -> GetPDG() <<"\t" << (*ScintStart1HC)[0] -> GetEkin()/MeV <<"\t" << (*ScintStart1HC)[0] -> GetEdep()/MeV << G4endl;
	//G4cout <<(*NaIHC)[0] -> GetPDG() <<"\t" << (*NaIHC)[0] -> GetEkin()/MeV <<"\t" << (*NaIHC)[0] -> GetEdep()/MeV << G4endl;


}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
