// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
#ifndef NTUPLEFILTER_TRACKISOLATIONFILTER_H
#define NTUPLEFILTER_TRACKISOLATIONFILTER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "Math/VectorUtil.h"

#include "TLorentzVector.h"
#include "TTree.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class decleration
//

class TrackIsolationFilter : public edm::EDFilter {
public:
     explicit TrackIsolationFilter (const edm::ParameterSet&);
     ~TrackIsolationFilter();

private:
  virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

  // ----------member data ---------------------------
  double dR_;
  double dzcut_;
  double dxycut_;
  double minPt_;

  double isoCut_;

  bool doTrkIsoVeto_;

  std::vector<int> exclPdgIdVec_;

  edm::InputTag pfCandidatesTag_;
  edm::InputTag vertexInputTag_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> PfcandTok_;
  edm::EDGetTokenT<edm::View<reco::Vertex> > VertexInputTok_;
  int vtxSize;
};

#endif

