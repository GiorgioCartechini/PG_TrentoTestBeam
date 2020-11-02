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

#include <stdlib.h>
#include "AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

AnalysisManager::AnalysisManager() 
{
}

AnalysisManager::~AnalysisManager() 
{
}

void AnalysisManager::CreateNtuples() 
{ 
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  G4cout << "Using " << manager->GetType() << G4endl;
  
  manager->SetVerboseLevel(1);
  manager->SetNtupleMerging(true);
  manager->SetFirstNtupleId(0);
  
  
  //creating a ntuple, containing the information about secondary particles
  manager -> CreateNtuple("Lyso", "Lyso");
  manager -> CreateNtupleIColumn("PDG");
  manager -> CreateNtupleIColumn("TargetVol");
  manager -> CreateNtupleDColumn("Ekin");
  //manager -> CreateNtupleDColumn("Time");
  manager -> FinishNtuple();

  manager -> CreateNtuple("NaI", "NaI");
  manager -> CreateNtupleIColumn("PDG");
  manager -> CreateNtupleIColumn("TargetVol");
  manager -> CreateNtupleDColumn("Ekin");
  //manager -> CreateNtupleDColumn("Time");
  manager -> FinishNtuple();

  
  manager -> CreateNtuple("ScintStart1", "ScintStart1");
  manager -> CreateNtupleIColumn("PDG");
  manager -> CreateNtupleIColumn("TargetVol");
  manager -> CreateNtupleDColumn("Ekin");
  //manager -> CreateNtupleDColumn("Time");
  manager -> FinishNtuple();

  manager -> CreateNtuple("ScintVeto1", "ScintVeto1");
  manager -> CreateNtupleIColumn("PDG");
  manager -> CreateNtupleIColumn("TargetVol");
  manager -> CreateNtupleDColumn("Ekin");
  //manager -> CreateNtupleDColumn("Time");
  manager -> FinishNtuple();

  manager -> CreateNtuple("ScintVeto2", "ScintVeto2");
  manager -> CreateNtupleIColumn("PDG");
  manager -> CreateNtupleIColumn("TargetVol");
  manager -> CreateNtupleDColumn("Ekin");
  //manager -> CreateNtupleDColumn("Time");
  manager -> FinishNtuple();

  manager -> CreateNtuple("ExitTarget", "ExitTarget");
  manager -> CreateNtupleIColumn("PDG");
  manager -> CreateNtupleIColumn("TargetVol");
  manager -> CreateNtupleDColumn("Ekin");
  //manager -> CreateNtupleDColumn("Time");
  manager -> FinishNtuple();

  manager -> CreateNtuple("LysoEdep", "LysoEdep");
  manager -> CreateNtupleDColumn("EDep");
  manager -> FinishNtuple();

  manager -> CreateNtuple("NaIEdep", "NaIEdep");
  manager -> CreateNtupleDColumn("EDep");
  manager -> FinishNtuple();


   
}

void AnalysisManager::OpenFile(G4String fileName)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  
  // Create a root file 
  //G4String fileName = "radioprotection";

  // Create directories  
  //manager->SetNtupleDirectoryName("radioprotection_ntuple");
  

  G4bool fileOpen = manager->OpenFile(fileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << fileName << G4endl;
    return;
  }
}

void AnalysisManager::FillLyso(G4int PDG, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleIColumn(0, 0, PDG);
  manager -> FillNtupleIColumn(0, 1, Volume);
  manager -> FillNtupleDColumn(0, 2, Ekin);
  //manager -> FillNtupleDColumn(0, 3, Edep);
  //manager -> FillNtupleDColumn(0, 3, Time);
  manager -> AddNtupleRow(0);  
}


void AnalysisManager::FillNaI(G4int PDG, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleIColumn(1, 0, PDG);
  manager -> FillNtupleIColumn(1, 1, Volume);
  manager -> FillNtupleDColumn(1, 2, Ekin);
  //manager -> FillNtupleDColumn(0, 3, Edep);
  //manager -> FillNtupleDColumn(0, 3, Time);
  manager -> AddNtupleRow(1);  
}

void AnalysisManager::FillScintStart1(G4int PDG, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleIColumn(2, 0, PDG);
  manager -> FillNtupleIColumn(2, 1, Volume);
  manager -> FillNtupleDColumn(2, 2, Ekin);
  //manager -> FillNtupleDColumn(0, 3, Edep);
  //manager -> FillNtupleDColumn(0, 3, Time);
  manager -> AddNtupleRow(2);  
}

void AnalysisManager::FillScintVeto1(G4int PDG, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleIColumn(3, 0, PDG);
  manager -> FillNtupleIColumn(3, 1, Volume);
  manager -> FillNtupleDColumn(3, 2, Ekin);
  //manager -> FillNtupleDColumn(0, 3, Edep);
  //manager -> FillNtupleDColumn(0, 3, Time);
  manager -> AddNtupleRow(3);  
}


void AnalysisManager::FillScintVeto2(G4int PDG, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleIColumn(4, 0, PDG);
  manager -> FillNtupleIColumn(4, 1, Volume);
  manager -> FillNtupleDColumn(4, 2, Ekin);
  //manager -> FillNtupleDColumn(0, 3, Edep);
  //manager -> FillNtupleDColumn(0, 3, Time);
  manager -> AddNtupleRow(4);  
}

void AnalysisManager::FillTarget(G4int PDG, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleIColumn(5, 0, PDG);
  manager -> FillNtupleIColumn(5, 1, Volume);
  manager -> FillNtupleDColumn(5, 2, Ekin);
  //manager -> FillNtupleDColumn(0, 3, Edep);
  //manager -> FillNtupleDColumn(0, 3, Time);
  manager -> AddNtupleRow(5);  
}

void AnalysisManager::FillLysoEdep(G4double Edep)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(6, 0, Edep);
  manager -> AddNtupleRow(6);  
}

void AnalysisManager::FillNaIEdep(G4double Edep)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(7, 0, Edep);
  manager -> AddNtupleRow(7);  
}
 
void AnalysisManager::finish() 
{   
    G4AnalysisManager* manager = G4AnalysisManager::Instance();    
    manager -> Write();
    manager -> CloseFile();  
}












