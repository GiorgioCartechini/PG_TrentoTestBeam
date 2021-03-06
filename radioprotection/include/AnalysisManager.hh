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


#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH 

#include "globals.hh"
#include "g4root.hh"

// Define the total number of columns in the ntuple
//const G4int MaxNtCol = 5;

class AnalysisManager
{ 

	public:
		AnalysisManager();
		~AnalysisManager();

		void CreateNtuples(); // booking the ROOT file

		void OpenFile(G4String fileName); //open output file

		void FillLyso(G4int PDG, G4double Time, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/);
		void FillNaI(G4int PDG, G4double Time, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/);
		void FillLaBr(G4int PDG, G4double Time, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/); 
		void FillScintStart1(G4int PDG, G4double Time, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/); 
		void FillScintVeto1(G4int PDG, G4double Time, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/); 
		void FillScintVeto2(G4int PDG, G4double Time, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/); 
		void FillTarget(G4int PDG, G4double Time, G4int Volume, G4double Ekin /*, G4double Edep, G4double Time*/);
		void FillLysoEdep(G4double Edep);
		void FillNaIEdep(G4double Edep); 
		void FillLaBrEdep(G4double Edep); 
		// Information about secondary particles

		void finish();
		// Close the ROOT file with all the results stored in nutples 

};

#endif




