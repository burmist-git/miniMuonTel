#ifndef RunAction_h
#define RunAction_h 1

//my
#include "HitDataStructure.hh"
#include "EventDataStructure.hh"

//G4
#include "G4Timer.hh"
#include "globals.hh"
#include "G4UserRunAction.hh"

//root
#include "TTree.h"

class TFile;
class G4Run;

class RunAction : public G4UserRunAction
{
public:

  RunAction();
  virtual ~RunAction();

public:
  virtual void BeginOfRunAction(const G4Run* aRun);
  virtual void EndOfRunAction(const G4Run* aRun);

public:
  void SetOutputFileName(G4String fileName) {_outputFileName = fileName;}
  G4String GetOutputFileName() { return _outputFileName;}
  TTree* tree;
  HitData HitInfo;
  EventData EventInfo;

  G4int _nPhot;
  static const G4int _nPhotMax = 2000000;
  G4int _TrackID[_nPhotMax];
  G4int _ParentID[_nPhotMax];
  G4double _Energy[_nPhotMax];
  G4double _Wavelength[_nPhotMax];
  G4double _Time[_nPhotMax];
  G4double _photPathLen[_nPhotMax];
  G4int _SecID[_nPhotMax];
  G4int _chID[_nPhotMax];
  G4double _PosX[_nPhotMax];
  G4double _PosY[_nPhotMax];
  G4double _PosZ[_nPhotMax];
  G4double _MomX[_nPhotMax];
  G4double _MomY[_nPhotMax];
  G4double _MomZ[_nPhotMax];
  G4double _trkMomX[_nPhotMax];
  G4double _trkMomY[_nPhotMax];
  G4double _trkMomZ[_nPhotMax];
  G4double _trkPosX[_nPhotMax];
  G4double _trkPosY[_nPhotMax];
  G4double _trkPosZ[_nPhotMax];
  G4double _trkT[_nPhotMax];
  G4double _trkLength[_nPhotMax];

private:
  G4Timer* timer;
  TFile* hfile;
  G4String _outputFileName;

};

#endif
