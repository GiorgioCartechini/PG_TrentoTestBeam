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

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "SteppingAction.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "EventAction.hh"

SteppingAction::SteppingAction(AnalysisManager* pAnalysis, EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)

{ 
analysis = pAnalysis;
fSecondary = 0; 
}

SteppingAction::~SteppingAction()
{ 
}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* Prestep = aStep->GetPreStepPoint();
  G4StepPoint* Poststep = aStep->GetPostStepPoint();

  G4String PhVolumeName =  Prestep->GetPhysicalVolume()-> GetName();

  //fPName =  aStep  ->  GetTrack() -> GetDefinition() -> GetParticleName();

    
    if(PhVolumeName == "physicalScintStart_1")
    {
      if(Prestep ->GetStepStatus() == fGeomBoundary)
      {
        fEventAction  -> SetEkinScStart1( Prestep -> GetKineticEnergy() );
        fEventAction  -> SetTimeScStart1( aStep->GetTrack() -> GetGlobalTime() );
        fEventAction  -> SetPDGScStart1( aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() );
        
        G4int vol = 0;
        if(aStep  ->  GetTrack() -> GetLogicalVolumeAtVertex()->GetName() == "logicTarget")
          vol = 1;
        else
          vol = 0;

        fEventAction  -> SetVertexVolScStart1(vol);

        analysis -> FillScintStart1(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode(), vol, Prestep -> GetKineticEnergy()/MeV);


        //fPName =  aStep  ->  GetTrack() -> GetDefinition() -> GetParticleName();
        //G4cout <<PhVolumeName <<"\t" << aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() << G4endl;
      }

      //if(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() == 22)
        //fEventAction  -> AddEdepScStart1(aStep -> GetTotalEnergyDeposit());
    }
    else if(PhVolumeName == "physicalScintVeto_1")
    {
      if(Prestep ->GetStepStatus() == fGeomBoundary)
      {
        fEventAction  -> SetEkinScVeto1( Prestep -> GetKineticEnergy() );
        fEventAction  -> SetTimeScVeto1( aStep->GetTrack() -> GetGlobalTime() );
        fEventAction  -> SetPDGScVeto1( aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() );

        G4int vol = 0;
        if(aStep  ->  GetTrack() -> GetLogicalVolumeAtVertex()->GetName() == "logicTarget")
          vol = 1;
        else
          vol = 0;

        fEventAction  -> SetVertexVolScVeto1(vol);

        analysis -> FillScintVeto1(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode(), vol, Prestep -> GetKineticEnergy()/MeV);

        //fPName =  aStep  ->  GetTrack() -> GetDefinition() -> GetParticleName();
        //G4cout <<volumeName <<"\t" << fPDG<<"\t" << fPName << "\t" << fTime/ns << G4endl;
        //G4cout <<PhVolumeName <<"\t" << vol <<"\t" << aStep  ->  GetTrack() -> GetLogicalVolumeAtVertex()->GetName() << G4endl;
      }

      //if(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() == 22)
        //fEventAction  -> AddEdepScVeto1(aStep -> GetTotalEnergyDeposit());
    }
    else if(PhVolumeName == "physicalScintVeto_2")
    {
      if(Prestep ->GetStepStatus() == fGeomBoundary)
      {
        fEventAction  -> SetEkinScVeto2( Prestep -> GetKineticEnergy() );
        fEventAction  -> SetTimeScVeto2( aStep->GetTrack() -> GetGlobalTime() );
        fEventAction  -> SetPDGScVeto2( aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() );

        G4int vol = 0;
        if(aStep  ->  GetTrack() -> GetLogicalVolumeAtVertex()->GetName() == "logicTarget")
          vol = 1;
        else
          vol = 0;

        fEventAction  -> SetVertexVolScVeto2(vol);

        analysis -> FillScintVeto2(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode(), vol, Prestep -> GetKineticEnergy()/MeV);

      }

      //if(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() == 22)
        //fEventAction  -> AddEdepScVeto2(aStep -> GetTotalEnergyDeposit());
    }
    else if(PhVolumeName == "physicalNaI")
    {
      if(Prestep ->GetStepStatus() == fGeomBoundary)
      {
        fEventAction  -> SetEkinNaI( Prestep -> GetKineticEnergy() );
        fEventAction  -> SetTimeNaI( aStep->GetTrack() -> GetGlobalTime() );
        fEventAction  -> SetPDGNaI( aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() );

        G4int vol = 0;
        if(aStep  ->  GetTrack() -> GetLogicalVolumeAtVertex()->GetName() == "logicTarget")
          vol = 1;
        else
          vol = 0;

        
        fEventAction  -> SetVertexVolNaI(vol);

        analysis -> FillNaI(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode(), vol, Prestep -> GetKineticEnergy()/MeV);

      }

      if(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() == 22)
        fEventAction  -> AddEdepNaI(aStep -> GetTotalEnergyDeposit());
    }
    else if(PhVolumeName == "physicalLyso")
    {
      if(Prestep ->GetStepStatus() == fGeomBoundary)
      {
        fEventAction  -> SetEkinLyso( Prestep -> GetKineticEnergy() );
        fEventAction  -> SetTimeLyso( aStep->GetTrack() -> GetGlobalTime() );
        fEventAction  -> SetPDGLyso( aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() );

        G4int vol = 0;
        if(aStep  ->  GetTrack() -> GetLogicalVolumeAtVertex()->GetName() == "logicTarget")
          vol = 1;
        else
          vol = 0;

        fEventAction  -> SetVertexVolLyso(vol);

        analysis -> FillLyso(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode(), vol, Prestep -> GetKineticEnergy()/MeV);
      }

        //fPName =  aStep  ->  GetTrack() -> GetDefinition() -> GetParticleName();
        //G4cout <<volumeName <<"\t" << fPDG<<"\t" << fPName << "\t" << fTime/ns << G4endl;
        if(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() == 22)
          fEventAction  -> AddEdepLyso(aStep -> GetTotalEnergyDeposit());
    }
      
    else if(Prestep->GetPhysicalVolume()-> GetName() == "physicalTarget")
    {
      if(Poststep ->GetStepStatus() == fGeomBoundary)
      {
        fEventAction  -> SetEkinTgt( Prestep -> GetKineticEnergy() );
        fEventAction  -> SetTimeTgt( aStep->GetTrack() -> GetGlobalTime() );
        fEventAction  -> SetPDGTgt( aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode() );

        G4int vol = 0;
        if(aStep  ->  GetTrack() -> GetLogicalVolumeAtVertex()->GetName() == "logicTarget")
          vol = 1;
        else
          vol = 0;

        fEventAction  -> SetVertexVolTgt(vol);

        analysis -> FillTarget(aStep->GetTrack() -> GetDynamicParticle() -> GetPDGcode(), vol, Poststep -> GetKineticEnergy()/MeV);

        //fPName =  aStep  ->  GetTrack() -> GetDefinition() -> GetParticleName();
        //G4cout <<volumeName <<"\t" << fPDG<<"\t" << fPName << "\t" << fTime/ns << G4endl;
      }
    }
}

