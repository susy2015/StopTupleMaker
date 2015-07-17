// =================================================================/uscms_data/d3/benwu/StopHad/CMSSW_7_3_0/src/SusyAnaTools/SkimsAUX/plugins/prodFatTops.cc
// 
//       Filename:  prodFatTops.cc
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  07/06/2015 11:22:24 AM
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu (benwu)
//          Email:  benwu@fnal.gov
//        Company:  Baylor University, CDF@FNAL, CMS@LPC
// 
// ===========================================================================

/**\class prodFatTops

Description: 

Implementation:
 
*/

#ifndef  __EDFilter_PRODFATTOPS_INC__
#define  __EDFilter_PRODFATTOPS_INC__

// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

// ROOT
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMS Top Tagger
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"


// HEP V2 Top Tagger
#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//
// class declaration
//

class prodFatTops : public edm::EDFilter {
  public:
    explicit prodFatTops(const edm::ParameterSet&);
    ~prodFatTops();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual bool beginRun(edm::Run&, edm::EventSetup const&);
    virtual bool endRun(edm::Run&, edm::EventSetup const&);
    virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

    std::vector<TLorentzVector> prodPatJet( edm::Handle<std::vector<pat::Jet> > &patjets) const;
    std::vector<double> prodAK4Btag() const;
    std::vector<TLorentzVector> prodPFJet( edm::Handle<std::vector<reco::PFJet> > &RecoJets) const;
    std::vector<TLorentzVector> prodGenJets(edm::Handle<std::vector<reco::GenJet> > &Genjets) const;

    std::vector<TLorentzVector> prodCMSTopTagger();
    std::vector<TLorentzVector> prodHEPV2TopTagger(edm::Handle<std::vector<reco::HTTTopJetTagInfo> > &HTTs) const;
    std::vector<TLorentzVector> prodSoftDropTagger(edm::Handle<std::vector<pat::Jet> > & SoftDropjets) const;
    std::vector<double> prodSoftDropTau32(edm::Handle<edm::ValueMap<float> > &SoftDroptau3, edm::Handle<edm::ValueMap<float> > &SoftDroptau2) const;
    std::vector<double> prodShowerDecontruction( edm::Handle<std::vector<reco::PFJet> > &RecoJets, edm::Handle<edm::ValueMap<double> > ShowerD_chi) const;
    // ----------member data ---------------------------

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Skiny jets ~~~~~
    edm::InputTag AK4Src_;
    edm::Handle<std::vector<pat::Jet> > AK4jets;

    std::string bTagKeyString_;

    edm::InputTag Gen4Src_;
    edm::Handle<std::vector<reco::GenJet> > Gen4jets;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fat Jets ~~~~~
    edm::InputTag Gen8Src_;
    edm::Handle<std::vector<reco::GenJet> > Gen8jets;

    edm::InputTag AK8Src_;
    edm::Handle<std::vector<pat::Jet> > AK8jets;

    edm::InputTag RecoAK8Src_;
    edm::Handle<std::vector<reco::PFJet> > RecoAK8jets;
    
    edm::InputTag HTTSrc_;
    edm::Handle<std::vector<reco::HTTTopJetTagInfo> > HTTjets;

    edm::InputTag ORHTTSrc_;
    edm::Handle<std::vector<reco::HTTTopJetTagInfo> > ORHTTjets;

    edm::InputTag AK8SoftDropSrc_;
    std::string AK8SoftDropTauSrc_;
    edm::Handle<std::vector<pat::Jet> > AK8SoftDropjets;
    edm::Handle<edm::ValueMap<float> > AK8SoftDropjets_tau2;
    edm::Handle<edm::ValueMap<float> > AK8SoftDropjets_tau3;

    edm::InputTag CA15SoftDropSrc_;
    std::string CA15SoftDropTauSrc_;
    edm::Handle<std::vector<pat::Jet> > CA15SoftDropjets;
    edm::Handle<edm::ValueMap<float> > CA15SoftDropjets_tau2;
    edm::Handle<edm::ValueMap<float> > CA15SoftDropjets_tau3;

    std::string ShowerD2Src_;
    edm::Handle<edm::ValueMap<double> > ShowerD_chi2;
    edm::Handle<edm::ValueMap<int> > ShowerD_nmj2;

    std::string ShowerD3Src_;
    edm::Handle<edm::ValueMap<double> > ShowerD_chi3;
    edm::Handle<edm::ValueMap<int> > ShowerD_nmj3;

    bool debug_;
};



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
prodFatTops::prodFatTops(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  // Skiny Jets
  AK4Src_         = iConfig.getParameter<edm::InputTag>("AK4Src");
  bTagKeyString_ = iConfig.getParameter<std::string>("bTagKeyString");
  Gen4Src_         = iConfig.getParameter<edm::InputTag>("Gen4Src");
  Gen8Src_         = iConfig.getParameter<edm::InputTag>("Gen8Src");

  // FAT Jets
  AK8Src_         = iConfig.getParameter<edm::InputTag>("AK8Src");
  RecoAK8Src_     = iConfig.getParameter<edm::InputTag>("RecoAK8Src");
  AK8SoftDropSrc_ = iConfig.getParameter<edm::InputTag>("ak8PFJetsSoftDrop");
  CA15SoftDropSrc_ = iConfig.getParameter<edm::InputTag>("ca15PFJetsSoftDrop");
  AK8SoftDropTauSrc_ = iConfig.getParameter<std::string>("ak8PFJetsSoftDropTau");
  CA15SoftDropTauSrc_ = iConfig.getParameter<std::string>("ca15PFJetsSoftDropTau");
  HTTSrc_         = iConfig.getParameter<edm::InputTag>("HTTSrc");
  ORHTTSrc_       = iConfig.getParameter<edm::InputTag>("ORHTTSrc");
  debug_          = iConfig.getParameter<bool>("debug");
  bTagKeyString_  = iConfig.getParameter<std::string>("bTagKeyString");
  ShowerD2Src_    = iConfig.getParameter<std::string>("ShowerD2Src");
  ShowerD3Src_    = iConfig.getParameter<std::string>("ShowerD3Src");

  produces<std::vector<TLorentzVector> >("jetsLVec");
  produces<std::vector<double> >("recoJetsBtag");

  produces<std::vector<double> >("AK8SD2Chi");
  produces<std::vector<double> >("AK8SD3Chi");

  produces<std::vector<double> >("AK8SoftDropTau32");
  produces<std::vector<double> >("CA15SoftDropTau32");
  
  produces<std::vector<TLorentzVector> >("Gen4LVec");
  produces<std::vector<TLorentzVector> >("Gen8LVec");
  produces<std::vector<TLorentzVector> >("CMSTopLVec");
  produces<std::vector<TLorentzVector> >("AK8LVec");
  produces<std::vector<TLorentzVector> >("HEPV2TopLVec");
  produces<std::vector<TLorentzVector> >("ORHEPV2TopLVec");
  produces<std::vector<TLorentzVector> >("AK8SoftDropLVec");
  produces<std::vector<TLorentzVector> >("CA15SoftDropLVec");

}


prodFatTops::~prodFatTops()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
  bool
prodFatTops::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  iEvent.getByLabel(AK4Src_, AK4jets);
  iEvent.getByLabel(Gen4Src_, Gen4jets);
  iEvent.getByLabel(Gen8Src_, Gen8jets);

  iEvent.getByLabel(AK8Src_, AK8jets);
  iEvent.getByLabel(RecoAK8Src_, RecoAK8jets);
  iEvent.getByLabel(HTTSrc_, HTTjets);
  iEvent.getByLabel(ORHTTSrc_, ORHTTjets);
  iEvent.getByLabel(AK8SoftDropSrc_, AK8SoftDropjets);
  iEvent.getByLabel(AK8SoftDropTauSrc_, "tau2", AK8SoftDropjets_tau2);
  iEvent.getByLabel(AK8SoftDropTauSrc_, "tau3", AK8SoftDropjets_tau3);
  iEvent.getByLabel(CA15SoftDropSrc_, CA15SoftDropjets);
  iEvent.getByLabel(CA15SoftDropTauSrc_, "tau2", CA15SoftDropjets_tau2);
  iEvent.getByLabel(CA15SoftDropTauSrc_, "tau3", CA15SoftDropjets_tau3);
  iEvent.getByLabel(ShowerD2Src_, "chi", ShowerD_chi2);
  iEvent.getByLabel(ShowerD3Src_, "chi", ShowerD_chi3);

  if (AK4jets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > AK4LVec(new std::vector<TLorentzVector>());
    std::auto_ptr<std::vector<double> > AK4Btag(new std::vector<double>());
    *AK4LVec.get() = prodPatJet(AK4jets);
    *AK4Btag.get() = prodAK4Btag();
    iEvent.put(AK4LVec, "jetsLVec");
    iEvent.put(AK4Btag, "recoJetsBtag");
  }


  if (Gen4jets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > Gen4LVec(new std::vector<TLorentzVector>());
    *Gen4LVec.get() = prodGenJets(Gen4jets);
    iEvent.put(Gen4LVec, "Gen4LVec" );
  }

  if (Gen8jets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > Gen8LVec(new std::vector<TLorentzVector>());
    *Gen8LVec.get() = prodGenJets(Gen8jets);
    iEvent.put(Gen8LVec, "Gen8LVec" );
  }


  if (RecoAK8jets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > AK8LVec(new std::vector<TLorentzVector>());
    *AK8LVec.get() = prodPFJet(RecoAK8jets);
    iEvent.put(AK8LVec, "AK8LVec");

    if (ShowerD_chi2.isValid())
    {
      std::auto_ptr<std::vector<double> > AK8SD2Chi(new std::vector<double>());
      *AK8SD2Chi.get() = prodShowerDecontruction(RecoAK8jets, ShowerD_chi2);
      iEvent.put(AK8SD2Chi, "AK8SD2Chi");
    }
    if (ShowerD_chi3.isValid())
    {
      std::auto_ptr<std::vector<double> > AK8SD3Chi(new std::vector<double>());
      *AK8SD3Chi.get() = prodShowerDecontruction(RecoAK8jets, ShowerD_chi3);
      iEvent.put(AK8SD3Chi, "AK8SD3Chi");
    }
  }

  if (AK8jets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > CMSTopLVec(new std::vector<TLorentzVector>());
    *CMSTopLVec.get() = prodCMSTopTagger();
    iEvent.put(CMSTopLVec, "CMSTopLVec");
  }


  if (HTTjets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > HEPV2TopLVec(new std::vector<TLorentzVector>());
    *HEPV2TopLVec.get() = prodHEPV2TopTagger(HTTjets);
    iEvent.put(HEPV2TopLVec, "HEPV2TopLVec");
  }

  if (ORHTTjets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > ORHEPV2TopLVec(new std::vector<TLorentzVector>());
    *ORHEPV2TopLVec.get() = prodHEPV2TopTagger(ORHTTjets);
    iEvent.put(ORHEPV2TopLVec, "ORHEPV2TopLVec");
  }

  if (AK8SoftDropjets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > SoftDropLVec(new std::vector<TLorentzVector>());
    *SoftDropLVec.get() = prodSoftDropTagger(AK8SoftDropjets);
    iEvent.put(SoftDropLVec, "AK8SoftDropLVec");
    if (AK8SoftDropjets_tau2.isValid() && AK8SoftDropjets_tau3.isValid())
    {
      std::auto_ptr<std::vector<double> > SoftDropTau32(new std::vector<double>());
      *SoftDropTau32.get() = prodSoftDropTau32(AK8SoftDropjets_tau3, AK8SoftDropjets_tau2);
      iEvent.put(SoftDropTau32, "AK8SoftDropTau32");
    }
  }

  if (CA15SoftDropjets.isValid())
  {
    std::auto_ptr<std::vector<TLorentzVector> > SoftDropLVec(new std::vector<TLorentzVector>());
    *SoftDropLVec.get() = prodSoftDropTagger(CA15SoftDropjets);
    iEvent.put(SoftDropLVec, "CA15SoftDropLVec");
    if (CA15SoftDropjets_tau2.isValid() && CA15SoftDropjets_tau3.isValid())
    {
      std::auto_ptr<std::vector<double> > SoftDropTau32(new std::vector<double>());
      *SoftDropTau32.get() = prodSoftDropTau32(CA15SoftDropjets_tau3, CA15SoftDropjets_tau2);
      iEvent.put(SoftDropTau32, "CA15SoftDropTau32");
    }
  }

  return true;
}

// ------------ method called once each job just before starting event loop  ------------
  void 
prodFatTops::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
prodFatTops::endJob() {
}

// ------------ method called when starting to processes a run  ------------
  bool 
prodFatTops::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
  bool 
prodFatTops::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
  bool 
prodFatTops::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
  bool 
prodFatTops::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void
prodFatTops::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodCMSTopTagger
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> prodFatTops::prodCMSTopTagger()
{
  std::vector<TLorentzVector> temp;
  for(unsigned int i=0; i < AK8jets->size(); ++i)
  {
    pat::Jet jet = AK8jets->at(i);

    reco::CATopJetTagInfo const * tagInfo =  dynamic_cast<reco::CATopJetTagInfo const *>( jet.tagInfo("caTop"));

    bool topTagged = false;
    if ( tagInfo != 0 ) {
      double minMass = tagInfo->properties().minMass;
      double topMass = tagInfo->properties().topMass;
      int nSubJets = tagInfo->properties().nSubJets;
      
      //std::cout << " min " << minMass <<" top " << topMass <<" nSub " << nSubJets << std::endl;

      if ( nSubJets > 2 && minMass > 50.0 && topMass > 150.0 )
        topTagged = true;
    }

    if (topTagged)
    {
      TLorentzVector perJetLVec;
      perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
      temp.push_back(perJetLVec);
    }


  }

  return temp;
}       // -----  end of function prodFatTops::prodCMSTopTagger  -----

// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodHEPV2TopTagger
//  Description:  /* cursor */
// ===========================================================================
std::vector<TLorentzVector> prodFatTops::prodHEPV2TopTagger(edm::Handle<std::vector<reco::HTTTopJetTagInfo> > &HTTs) const
{
  std::vector<TLorentzVector> temp;
  for(unsigned int i=0; i < HTTs->size(); ++i)
  {
    reco::HTTTopJetTagInfo HTTInfo = HTTs->at(i);
    reco::HTTTopJetProperties HTTProp = HTTInfo.properties();
    TLorentzVector perJetLVec;
    perJetLVec.SetPtEtaPhiM(HTTProp.fjPt, HTTProp.fjEta, HTTProp.fjPhi, HTTProp.fjMass);
    temp.push_back(perJetLVec);
  }
  return temp;
}       // -----  end of function prodFatTops::prodHEPV2TopTagger  -----

// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodSoftDropTagger
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> prodFatTops::prodSoftDropTagger(edm::Handle<std::vector<pat::Jet> > & SoftDropjets) const
{
  std::vector<TLorentzVector> temp;
  for(unsigned int i=0; i < SoftDropjets->size(); ++i)
  {
    pat::Jet jet = SoftDropjets->at(i);

    //double tau2 = jet.userFloat("NjettinessAK8CHS:tau2");    //  Access the n-subjettiness variables
    //double tau3 = jet.userFloat("NjettinessAK8CHS:tau3");    // 

    //double tau32 = tau3/tau2; 
    //std::cout <<  " tau 2 " << tau2 <<" tau3 " << tau3 << " tau32 " << tau32 << std::endl;

    if (jet.mass() > 110 && jet.mass() < 210 )
    //if (jet.mass() > 110 && jet.mass() < 210 && tau32 < 0.6)
    {
      TLorentzVector perJetLVec;
      perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
      temp.push_back(perJetLVec);
    }
  }

  return temp;
}       // -----  end of function prodFatTops::prodSoftDropTagger  -----

// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodSoftDropTau32
//  Description:  /* cursor */
// ===========================================================================
std::vector<double> prodFatTops::prodSoftDropTau32(edm::Handle<edm::ValueMap<float> > &SoftDroptau3,
    edm::Handle<edm::ValueMap<float> > &SoftDroptau2) const
{
  std::vector<double> temp;
  for(unsigned int i=0; i < SoftDroptau2->size(); ++i)
  {
    double out = -999;
    if (SoftDroptau2->get(i) != 0)
    {
     out =  SoftDroptau3->get(i) / SoftDroptau2->get(i);
    }

    temp.push_back(out);
  }

  return temp;
}       // -----  end of function prodFatTops::prodSoftDropTau32  -----

// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodAK4CHS
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> prodFatTops::prodPatJet( edm::Handle<std::vector<pat::Jet> > &patjets) const
{
  
  std::vector<TLorentzVector> temp;
  for(unsigned int i=0; i < patjets->size(); ++i)
  {
    pat::Jet jet = patjets->at(i);
    if (jet.pt() > 20)
    {
      TLorentzVector perJetLVec;
      perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
      temp.push_back(perJetLVec);
    }
  }

  return temp;
}       // -----  end of function prodFatTops::prodAK4CHS  -----

// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodAK4Btag
//  Description:  /* cursor */
// ===========================================================================
std::vector<double> prodFatTops::prodAK4Btag() const
{
  std::vector<double> temp;
  for(unsigned int i=0; i < AK4jets->size(); ++i)
  {
    pat::Jet jet = AK4jets->at(i);
    if (jet.pt() > 20)
    {
      temp.push_back(jet.bDiscriminator(bTagKeyString_.c_str()));
    }
  }

  return temp;
}       // -----  end of function prodFatTops::prodAK4Btag  -----

// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodGenJets
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> prodFatTops::prodGenJets(edm::Handle<std::vector<reco::GenJet> > &Genjets) const
{
  std::vector<TLorentzVector> temp;
  for(unsigned int i=0; i < Genjets->size(); ++i)
  {
    reco::GenJet jet = Genjets->at(i);
    if (jet.pt() > 10)
    {
      TLorentzVector perJetLVec;
      perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
      temp.push_back(perJetLVec);
    }
  }

  return temp;
}       // -----  end of function prodFatTops::prodGenJets  -----

// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodPFJet
//  Description:  /* cursor */
// ===========================================================================
std::vector<TLorentzVector> prodFatTops::prodPFJet( edm::Handle<std::vector<reco::PFJet> > &RecoJets) const
{
  std::vector<TLorentzVector> temp;
  for(unsigned int i=0; i < RecoJets->size(); ++i)
  {
    reco::PFJet jet = RecoJets->at(i);
      TLorentzVector perJetLVec;
      perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
      temp.push_back(perJetLVec);
  }
  return temp;
}       // -----  end of function prodFatTops::prodPFJet  -----


// ===  FUNCTION  ============================================================
//         Name:  prodFatTops::prodShowerDecontruction
//  Description:  [> cursor <]
// ===========================================================================
std::vector<double> prodFatTops::prodShowerDecontruction( edm::Handle<std::vector<reco::PFJet> > &RecoJets, edm::Handle<edm::ValueMap<double> > ShowerD_chi) const
{
  std::vector<double> temp;
  assert(RecoJets->size() == ShowerD_chi->size());

  for(unsigned int i=0; i < RecoJets->size(); ++i)
  {
    double out = -999;
    if (ShowerD_chi->get(i) != 0)
    {
      out = log(ShowerD_chi->get(i));
    }
    temp.push_back(out);
  }
  return temp;
}       // -----  end of function prodFatTops::prodShowerDecontruction  -----

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(prodFatTops);
#endif   // ----- #ifndef __EDFilter_PRODFATTOPS_INC__  -----



