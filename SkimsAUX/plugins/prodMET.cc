#include <memory>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"

#include "TLorentzVector.h"

class prodMET : public edm::EDFilter {

  public:

    explicit prodMET(const edm::ParameterSet & iConfig);
    ~prodMET();

  private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

    edm::InputTag metSrc_;
    bool debug_;
    bool addcalomet_;

    bool isData_;

    edm::EDGetTokenT<edm::View<pat::MET> >MetTok_;

    std::vector<pat::MET::METUncertainty> uncUpList, uncDownList;
};


prodMET::prodMET(const edm::ParameterSet & iConfig) {

  isData_ = true;

  metSrc_      = iConfig.getParameter<edm::InputTag>("metSrc");

  debug_       = iConfig.getParameter<bool>("debug");
  addcalomet_  = iConfig.getParameter<bool>("addcalomet");

  MetTok_  = consumes<edm::View<pat::MET> >(metSrc_);

  uncUpList = {pat::MET::JetResUp, pat::MET::JetEnUp, pat::MET::MuonEnUp, pat::MET::ElectronEnUp, pat::MET::TauEnUp, pat::MET::UnclusteredEnUp, pat::MET::PhotonEnUp};
  uncDownList = {pat::MET::JetResDown, pat::MET::JetEnDown, pat::MET::MuonEnDown, pat::MET::ElectronEnDown, pat::MET::TauEnDown, pat::MET::UnclusteredEnDown, pat::MET::PhotonEnDown};

  produces<float>("met");
  produces<float>("metphi");
  produces<float>("genmet");
  produces<float>("genmetphi");
  produces<float>("calomet");
  produces<float>("calometphi");
  produces<std::vector<float> >("metMagUp");
  produces<std::vector<float> >("metMagDown");
  produces<std::vector<float> >("metPhiUp");
  produces<std::vector<float> >("metPhiDown");
}


prodMET::~prodMET() {
}


bool prodMET::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if( !iEvent.isRealData() ) isData_ = false;

  edm::Handle<edm::View<pat::MET> > met;
  iEvent.getByToken(MetTok_, met);

  std::unique_ptr<float> metPtr(new float);
  std::unique_ptr<float> metphiPtr(new float);

  std::unique_ptr<float> genmetPtr(new float);
  std::unique_ptr<float> genmetphiPtr(new float);

  std::unique_ptr<float> calometPtr(new float);
  std::unique_ptr<float> calometphiPtr(new float);

  std::unique_ptr<std::vector<float> > metMagUp_ (new std::vector<float>(uncUpList.size(), 0.));
  std::unique_ptr<std::vector<float> > metPhiUp_ (new std::vector<float>(uncUpList.size(), 0.));
  std::unique_ptr<std::vector<float> > metMagDown_ (new std::vector<float>(uncDownList.size(), 0.));
  std::unique_ptr<std::vector<float> > metPhiDown_ (new std::vector<float>(uncDownList.size(), 0.));

  *metPtr = (*met)[0].pt();
  *metphiPtr = (*met)[0].phi();

  for(unsigned int iu=0; iu<uncUpList.size(); ++iu){
     (*metMagUp_)[iu] = met->at(0).shiftedPt(uncUpList[iu], pat::MET::Type1);
     (*metPhiUp_)[iu] = met->at(0).shiftedPhi(uncUpList[iu], pat::MET::Type1);
  }
  for(unsigned int iu=0; iu<uncDownList.size(); ++iu){
     (*metMagDown_)[iu] = met->at(0).shiftedPt(uncDownList[iu], pat::MET::Type1);
     (*metPhiDown_)[iu] = met->at(0).shiftedPhi(uncDownList[iu], pat::MET::Type1);
  }

  if( met.isValid() && !isData_ ){
     if( met->at(0).genMET() ){
        const reco::GenMET * theGenMET(met->at(0).genMET());
        *genmetPtr = theGenMET->pt();
        *genmetphiPtr = theGenMET->phi();
     }
  }

  if( met.isValid() && addcalomet_ ){
     *calometPtr = met->at(0).caloMETPt();
     *calometphiPtr = met->at(0).caloMETPhi();
  }

  iEvent.put(std::move(metMagUp_), "metMagUp");
  iEvent.put(std::move(metPhiUp_), "metPhiUp");
  iEvent.put(std::move(metMagDown_), "metMagDown");
  iEvent.put(std::move(metPhiDown_), "metPhiDown");
  iEvent.put(std::move(metPtr), "met");
  iEvent.put(std::move(metphiPtr), "metphi");
  iEvent.put(std::move(genmetPtr), "genmet");
  iEvent.put(std::move(genmetphiPtr), "genmetphi");
  if( addcalomet_ ){
    iEvent.put(std::move(calometPtr), "calomet");
    iEvent.put(std::move(calometphiPtr), "calometphi");
  }

  return true;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(prodMET);
