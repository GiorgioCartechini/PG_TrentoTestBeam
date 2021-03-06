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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "SensitiveDetectorHit.hh"

#include "AnalysisManager.hh"

#include "globals.hh"

#include "G4Event.hh"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

class EventAction : public G4UserEventAction
{
public:
  EventAction(AnalysisManager* analysis);
  virtual ~EventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);

  void SetPDGScStart1     (G4int pdg)      { fPDGScStart1 = pdg; };
  void SetVertexVolScStart1     (G4int volume)      { fVolScStart1 = volume; };
  void SetEkinScStart1     (G4double EndEkin)      { fEkinScStart1 = EndEkin; };
  void AddEdepScStart1     (G4double de)      { fEdepScStart1 += de; };
  void SetTimeScStart1     (G4double t)      { fTimeScStart1 = t; };

  void SetPDGScVeto1     (G4int pdg)      { fPDGScVeto1 = pdg; };
  void SetVertexVolScVeto1     (G4int volume)      { fVolScVeto1 = volume; };
  void SetEkinScVeto1     (G4double EndEkin)      { fEkinScVeto1 = EndEkin; };
  void AddEdepScVeto1     (G4double de)      { fEdepScVeto1 += de; };
  void SetTimeScVeto1     (G4double t)      { fTimeScVeto1 = t; };

  void SetPDGScVeto2     (G4int pdg)      { fPDGScVeto2 = pdg; };
  void SetVertexVolScVeto2     (G4int volume)      { fVolScVeto2 = volume; };
  void SetEkinScVeto2     (G4double EndEkin)      { fEkinScVeto2 = EndEkin; };
  void AddEdepScVeto2     (G4double de)      { fEdepScVeto2 += de; };
  void SetTimeScVeto2     (G4double t)      { fTimeScVeto2 = t; };

  void SetPDGNaI     (G4int pdg)      { fPDGNaI = pdg; };
  void SetVertexVolNaI     (G4int volume)      { fVolNaI = volume; };
  void SetEkinNaI     (G4double EndEkin)      { fEkinNaI = EndEkin; };
  void AddEdepNaI     (G4double de)      { fEdepNaI += de; };
  void SetTimeNaI     (G4double t)      { fTimeNaI = t; };

  void SetPDGLyso     (G4int pdg)      { fPDGLyso = pdg; };
  void SetVertexVolLyso    (G4int volume)      { fVolLyso = volume; };
  void SetEkinLyso     (G4double EndEkin)      { fEkinLyso = EndEkin; };
  void AddEdepLyso     (G4double de)      { fEdepLyso += de; };
  void SetTimeLyso     (G4double t)      { fTimeLyso = t; };

  void SetPDGLaBr     (G4int pdg)      { fPDGLaBr = pdg; };
  void SetVertexVolLaBr    (G4int volume)      { fVolLaBr = volume; };
  void SetEkinLaBr     (G4double EndEkin)      { fEkinLaBr = EndEkin; };
  void AddEdepLaBr     (G4double de)      { fEdepLaBr += de; };
  void SetTimeLaBr     (G4double t)      { fTimeLaBr = t; };


  void SetPDGTgt     (G4int pdg)      { fPDGTgt = pdg; };
  void SetVertexVolTgt    (G4int volume)      { fVolTgt = volume; };
  void SetEkinTgt     (G4double EndEkin)      { fEkinTgt = EndEkin; };
  void AddEdepTgt     (G4double de)      { fEdepTgt += de; };
  void SetTimeTgt     (G4double t)      { fTimeTgt = t; };
    
    
  //Get Methods
  G4int GetPDGScStart1() const      { return fPDGScStart1; };
  G4int GetVertexVolScStart1() const      { return fVolScStart1; };
  G4double GetEkinScStart1() const     { return fEkinScStart1; };
  G4double GetEdepScStart1() const     { return fEdepScStart1; };
  G4double GetTimeScStart1() const     { return fTimeScStart1; };
  
  G4int GetPDGScVeto1() const      { return fPDGScVeto1; };
  G4int GetVertexVolScVeto1() const      { return fVolScVeto1; };
  G4double GetEkinScVeto1() const     { return fEkinScVeto1; };
  G4double GetEdepScVeto1() const     { return fEdepScVeto1; };
  G4double GetTimeScVeto1() const     { return fTimeScVeto1; };

  G4int GetPDGScVeto2() const      { return fPDGScVeto2; };
  G4int GetVertexVolScVeto2() const      { return fVolScVeto2; };
  G4double GetEkinScVeto2() const     { return fEkinScVeto2; };
  G4double GetEdepScVeto2() const     { return fEdepScVeto2; };
  G4double GetTimeScVeto2() const     { return fTimeScVeto2; };

  G4int GetPDGNaI() const      { return fPDGNaI; };
  G4int GetVertexVolNaI() const      { return fVolNaI; };
  G4double GetEkinNaI() const     { return fEkinNaI; };
  G4double GetEdepNaI() const     { return fEdepNaI; };
  G4double GetTimeNaI() const     { return fTimeNaI; };

  G4int GetPDGLyso() const      { return fPDGLyso; };
  G4int GetVertexVolLyso() const      { return fVolLyso; };
  G4double GetEkinLyso() const     { return fEkinLyso; };
  G4double GetEdepLyso() const     { return fEdepLyso; };
  G4double GetTimeLyso() const     { return fTimeLyso; };

  G4int GetPDGLaBr() const      { return fPDGLaBr; };
  G4int GetVertexVolLaBr() const      { return fVolLaBr; };
  G4double GetEkinLaBr() const     { return fEkinLaBr; };
  G4double GetEdepLaBr() const     { return fEdepLaBr; };
  G4double GetTimeLaBr() const     { return fTimeLaBr; };


  G4int GetPDGTgt() const      { return fPDGTgt; };
  G4int GetVertexVolTgt() const      { return fVolTgt; };
  G4double GetEkinTgt() const     { return fEkinTgt; };
  G4double GetEdepTgt() const     { return fEdepTgt; };
  G4double GetTimeTgt() const     { return fTimeTgt; };


private:
  // methods
  SensitiveDetectorHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
                            
  AnalysisManager* analysisMan;
  
  // data members                   
  G4int  fLysoHCID;
  G4int  fNaIHCID;
  G4int  fScintStart1HCID;
  G4int  fScintVeto1HCID;
  G4int  fScintVeto2HCID;


  G4int     fPDGScStart1;
  G4int     fVolScStart1;
  G4double  fEkinScStart1;
  G4double  fEdepScStart1;
  G4double  fTimeScStart1;

  G4int     fPDGScVeto1;
  G4int     fVolScVeto1;
  G4double  fEkinScVeto1;
  G4double  fEdepScVeto1;
  G4double  fTimeScVeto1;

  G4int     fPDGScVeto2;
  G4int     fVolScVeto2;
  G4double  fEkinScVeto2;
  G4double  fEdepScVeto2;
  G4double  fTimeScVeto2;

  G4int     fPDGNaI;
  G4int     fVolNaI;
  G4double  fEkinNaI;
  G4double  fEdepNaI;
  G4double  fTimeNaI;

  G4int     fPDGLyso;
  G4int     fVolLyso;
  G4double  fEkinLyso;
  G4double  fEdepLyso;
  G4double  fTimeLyso;

  G4int     fPDGLaBr;
  G4int     fVolLaBr;
  G4double  fEkinLaBr;
  G4double  fEdepLaBr;
  G4double  fTimeLaBr;


  G4int     fPDGTgt;
  G4int     fVolTgt;
  G4double  fEkinTgt;
  G4double  fEdepTgt;
  G4double  fTimeTgt;
};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
