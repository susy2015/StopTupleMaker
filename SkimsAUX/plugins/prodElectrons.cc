#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
//#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "StopTupleMaker/Skims/plugins/ElectronEffectiveArea.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/METReco/interface/MET.h"

#include "TLorentzVector.h"

#include "StopTupleMaker/SkimsAUX/plugins/common.h" 

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;

class prodElectrons : public edm::EDFilter 
{
  enum elesIDLevel {VETO, LOOSE, MEDIUM, TIGHT};
 public:
   explicit prodElectrons(const edm::ParameterSet & iConfig);
   ~prodElectrons();
 private:
  virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

  edm::InputTag electronSrc_;
  edm::InputTag conversionsSrc_;
  edm::InputTag vtxSrc_;
  edm::InputTag metSrc_;
  edm::InputTag beamSpotSrc_;
  edm::InputTag pfCandsSrc_;
  edm::InputTag rhoSrc_;
  edm::EDGetTokenT< edm::View<pat::Electron> > ElecTok_;
  edm::EDGetTokenT< std::vector<reco::Vertex> > VtxTok_;
  edm::EDGetTokenT<edm::View<reco::MET> > MetTok_;
  edm::EDGetTokenT< std::vector<reco::Conversion> > ConversionsTok_;
  edm::EDGetTokenT<reco::BeamSpot> BeamSpotTok_;
  edm::EDGetTokenT<pat::PackedCandidateCollection>  PfcandTok_;
  edm::EDGetTokenT<double> RhoTok_;

  edm::InputTag vetoElectronID;
  edm::InputTag looseElectronID;
  edm::InputTag mediumElectronID;
  edm::InputTag tightElectronID;

  edm::InputTag vetoID;
  edm::InputTag looseID;
  edm::InputTag mediumID;
  edm::InputTag tightID;

  edm::EDGetTokenT<edm::ValueMap<bool> >  vetoElectronIDTok_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  looseElectronIDTok_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  mediumElectronIDTok_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  tightElectronIDTok_;

  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > vetoIdFullInfoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > looseIdFullInfoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > mediumIdFullInfoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > tightIdFullInfoMapToken_;

  bool doEleVeto_, dod0dz_;
  int doEleIso_; // 0: don't do any isolation; 1: relIso;  2: miniIso
  double minElePt_, maxEleEta_;
  bool debug_;
  double minElePtForElectron2Clean_, maxEleMiniIso_;

  bool passElectronID(const pat::Electron & ele, const edm::Handle< std::vector<reco::Vertex> > & vertices, const elesIDLevel level);
  bool passElectronISO(const pat::Electron & ele, const double relIso, const elesIDLevel level);
};


typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;


prodElectrons::prodElectrons(const edm::ParameterSet & iConfig)
{
  electronSrc_   = iConfig.getParameter<edm::InputTag>("ElectronSource");
  conversionsSrc_= iConfig.getParameter<edm::InputTag>("ConversionsSource");
  vtxSrc_        = iConfig.getParameter<edm::InputTag>("VertexSource");
  metSrc_        = iConfig.getParameter<edm::InputTag>("metSource");
  beamSpotSrc_   = iConfig.getParameter<edm::InputTag>("BeamSpotSource");
  pfCandsSrc_    = iConfig.getParameter<edm::InputTag>("PFCandSource");
  minElePt_      = iConfig.getParameter<double>("MinElePt");
  maxEleEta_     = iConfig.getParameter<double>("MaxEleEta");
  doEleVeto_     = iConfig.getParameter<bool>("DoElectronVeto");
  dod0dz_        = iConfig.getParameter<bool>("Dod0dz");
  doEleIso_      = iConfig.getParameter<int>("DoElectronIsolation");
  maxEleMiniIso_ = iConfig.getParameter<double>("MaxEleMiniIso");
  debug_         = iConfig.getParameter<bool>("Debug");
  rhoSrc_       = iConfig.getParameter<edm::InputTag>("RhoSource");

  vetoElectronID = iConfig.getParameter<edm::InputTag>("VetoElectronID");
  looseElectronID = iConfig.getParameter<edm::InputTag>("LooseElectronID");
  mediumElectronID = iConfig.getParameter<edm::InputTag>("MediumElectronID");
  tightElectronID = iConfig.getParameter<edm::InputTag>("TightElectronID");

  vetoID = iConfig.getParameter<edm::InputTag>("VetoElectronID");
  looseID = iConfig.getParameter<edm::InputTag>("LooseElectronID");
  mediumID = iConfig.getParameter<edm::InputTag>("MediumElectronID");
  tightID = iConfig.getParameter<edm::InputTag>("TightElectronID");

  minElePtForElectron2Clean_ = iConfig.getUntrackedParameter<double>("minElePtForElectron2Clean", 10);
  
  ConversionsTok_ = consumes< std::vector<reco::Conversion> >(conversionsSrc_);
  BeamSpotTok_ = consumes<reco::BeamSpot>(beamSpotSrc_);
  ElecTok_ = consumes< edm::View<pat::Electron> >(electronSrc_);
  MetTok_ = consumes<edm::View<reco::MET>>(metSrc_);
  VtxTok_ = consumes<std::vector<reco::Vertex>>(vtxSrc_);
  PfcandTok_ = consumes<pat::PackedCandidateCollection>(pfCandsSrc_);
  RhoTok_ = consumes<double>(rhoSrc_);

  vetoElectronIDTok_=consumes<edm::ValueMap<bool> >(vetoElectronID);
  looseElectronIDTok_=consumes<edm::ValueMap<bool> >(looseElectronID);
  mediumElectronIDTok_=consumes<edm::ValueMap<bool> >(mediumElectronID);
  tightElectronIDTok_=consumes<edm::ValueMap<bool> >(tightElectronID);

  vetoIdFullInfoMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >  (vetoID);
  looseIdFullInfoMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> > (looseID);
  mediumIdFullInfoMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> > (mediumID);
  tightIdFullInfoMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >  (tightID);

  produces<std::vector<pat::Electron> >("");
  produces<std::vector<pat::Electron> >("ele2Clean");  

  produces<std::vector<int> >("elesFlagVeto");
  produces<std::vector<int> >("elesFlagMedium");
  produces<std::vector<int> >("elesFlagLoose");
  produces<std::vector<int> >("elesFlagTight");
  produces<std::vector<TLorentzVector> >("elesLVec");
  produces<std::vector<float> >("elesCharge");
  produces<std::vector<float> >("elesMtw");
  produces<std::vector<float> >("elesRelIso");
  produces<std::vector<bool> >("elesisEB");
  produces<std::vector<float> >("elesMiniIso");
  produces<std::vector<float> >("elespfActivity");
  produces<int>("nElectrons");
 
  produces< std::vector< bool > >("vetoElectronID");
  produces< std::vector< bool > >("looseElectronID");
  produces< std::vector< bool > >("mediumElectronID");
  produces< std::vector< bool > >("tightElectronID");

  produces< std::vector< bool > >("vetoID");
  produces< std::vector< bool > >("looseID");
  produces< std::vector< bool > >("mediumID");
  produces< std::vector< bool > >("tightID");

}


prodElectrons::~prodElectrons()
{
}


bool prodElectrons::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  // electrons
  edm::Handle< edm::View<pat::Electron> > electrons;   
  iEvent.getByToken(ElecTok_, electrons);

  // conversions
  edm::Handle< std::vector<reco::Conversion> > conversions;
  iEvent.getByToken(ConversionsTok_, conversions);

  // beam spot
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(BeamSpotTok_, beamspot);
  //const reco::BeamSpot &beamSpot = *(beamspot.product());

  edm::Handle<edm::ValueMap<bool> >   veto_id_decisions_;
  iEvent.getByToken(vetoElectronIDTok_,veto_id_decisions_);
  edm::Handle<edm::ValueMap<bool> >   loose_id_decisions_;
  iEvent.getByToken(looseElectronIDTok_,loose_id_decisions_);
  edm::Handle<edm::ValueMap<bool> >   medium_id_decisions_;
  iEvent.getByToken(mediumElectronIDTok_,medium_id_decisions_);
  edm::Handle<edm::ValueMap<bool> >   tight_id_decisions_;
  iEvent.getByToken(tightElectronIDTok_,tight_id_decisions_); 
 
  // vertices
  edm::Handle< std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(VtxTok_, vertices);

  edm::Handle<edm::View<reco::MET> > met;
  iEvent.getByToken(MetTok_, met);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(PfcandTok_, pfcands);

  edm::Handle< double > rho_;
  //iEvent.getByLabel("fixedGridRhoFastjetCentralNeutral", rho_); // Central rho recommended for SUSY
  iEvent.getByToken(RhoTok_,rho_);
  double rho = *rho_;

  // check which ones to keep
  std::unique_ptr<std::vector<pat::Electron> > prod(new std::vector<pat::Electron>());
  std::unique_ptr<std::vector<pat::Electron> > ele2Clean(new std::vector<pat::Electron>());

  std::unique_ptr<std::vector<TLorentzVector> > elesLVec(new std::vector<TLorentzVector>());
  std::unique_ptr<std::vector<float> > elesCharge(new std::vector<float>());
  std::unique_ptr<std::vector<float> > elesMtw(new std::vector<float>());
  std::unique_ptr<std::vector<float> > elesRelIso(new std::vector<float>());
  std::unique_ptr<std::vector<bool> > elesisEB(new std::vector<bool>());
  std::unique_ptr<std::vector<float> > elesMiniIso(new std::vector<float>());
  std::unique_ptr<std::vector<float> > elespfActivity(new std::vector<float>());

  std::unique_ptr<std::vector<int> > elesFlagVeto(new std::vector<int>());
  std::unique_ptr<std::vector<int> > elesFlagMedium(new std::vector<int>());
  std::unique_ptr<std::vector<int> > elesFlagLoose(new std::vector<int>());
  std::unique_ptr<std::vector<int> > elesFlagTight(new std::vector<int>());

  edm::Handle<edm::ValueMap<vid::CutFlowResult> > veto_id_cutflow_;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > loose_id_cutflow_;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > tight_id_cutflow_;

  iEvent.getByToken(vetoIdFullInfoMapToken_,veto_id_cutflow_);
  iEvent.getByToken(looseIdFullInfoMapToken_,loose_id_cutflow_);
  iEvent.getByToken(mediumIdFullInfoMapToken_,medium_id_cutflow_);
  iEvent.getByToken(tightIdFullInfoMapToken_,tight_id_cutflow_);

  auto   Electron_vetoID  = std::make_unique<std::vector<bool>>();
  auto   Electron_looseID  = std::make_unique<std::vector<bool>>();
  auto   Electron_mediumID  = std::make_unique<std::vector<bool>>();
  auto   Electron_tightID  = std::make_unique<std::vector<bool>>();

  // loop on electrons
  for( edm::View<pat::Electron>::const_iterator ele = electrons->begin(); ele != electrons->end(); ele++ )
  {

    const edm::Ptr<pat::Electron> elePtr(electrons, ele - electrons->begin() );
    bool passveto = (*veto_id_decisions_)[ elePtr ];
    bool passloose = (*loose_id_decisions_)[ elePtr ];
    bool passmedium = (*medium_id_decisions_)[ elePtr ];
    bool passtight = (*tight_id_decisions_)[ elePtr ];

    vid::CutFlowResult mediumIdIsoMasked = (*medium_id_cutflow_)[ elePtr ].getCutFlowResultMasking("GsfEleEffAreaPFIsoCut_0");
    bool iPassMediumIDOnly_ = mediumIdIsoMasked.cutFlowPassed();
    
    vid::CutFlowResult looseIdIsoMasked = (*loose_id_cutflow_)[ elePtr ].getCutFlowResultMasking("GsfEleEffAreaPFIsoCut_0");
    bool iPassLooseIDOnly_ = looseIdIsoMasked.cutFlowPassed();
    
    vid::CutFlowResult vetoIdIsoMasked = (*veto_id_cutflow_)[ elePtr ].getCutFlowResultMasking("GsfEleEffAreaPFIsoCut_0");
    bool iPassVetoIDOnly_ = vetoIdIsoMasked.cutFlowPassed();
    
    vid::CutFlowResult tightIdIsoMasked = (*tight_id_cutflow_)[ elePtr ].getCutFlowResultMasking("GsfEleEffAreaPFIsoCut_0");
    bool iPassTightIDOnly_ = tightIdIsoMasked.cutFlowPassed();

    float pt = ele->pt();
    if (ele->pt() < minElePt_) continue;

    // get the ID variables from the electron object
    // kinematic variables
    bool isEB = ele->isEB() ? true : false;

    bool isVetoID = passElectronID((*ele), vertices, VETO);
    bool isLooseID = passElectronID((*ele), vertices, LOOSE);
    bool isMediumID = passElectronID((*ele), vertices, MEDIUM);
    bool isTightID = passElectronID((*ele), vertices, TIGHT);

    if( ! (isVetoID || isMediumID) ) continue;

    // isolation cuts                                                                                                                                        
    reco::GsfElectron::PflowIsolationVariables pfIso = ele->pfIsolationVariables();
    float absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );

    // compute final isolation
    float iso = absiso/pt;
    //double miniIso = commonFunctions::getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*ele)), 0.05, 0.2, 10., false, false);
    float miniIso = commonFunctions::GetMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*ele)), "electron", rho);
    float pfActivity = commonFunctions::GetMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*ele)), "electron", rho, true);

    if(doEleIso_ == 1 ) 
    {
      if( !passElectronISO((*ele), iso, VETO) ) continue;
    } 
    else if(doEleIso_ == 2 )
    {
      if( miniIso >= maxEleMiniIso_ ) continue;
    }

    // electron is ID'd and isolated! - only accept if vertex present
    if (vertices->size()>0)
    {
      prod->push_back(pat::Electron(*ele));
      TLorentzVector perLVec; perLVec.SetPtEtaPhiE(ele->pt(), ele->eta(), ele->phi(), ele->energy());
      elesLVec->push_back(perLVec);

      float mtw = sqrt( 2*( (*met)[0].pt()*ele->pt() -( (*met)[0].px()*ele->px() + (*met)[0].py()*ele->py() ) ) );

      elesCharge->push_back(ele->charge());
      elesMtw->push_back(mtw);
      elesRelIso->push_back(iso);
      elesisEB->push_back(isEB);
      elesMiniIso->push_back(miniIso);

      if( isVetoID ) elesFlagVeto->push_back(1); else elesFlagVeto->push_back(0); 
      if( isMediumID ) elesFlagMedium->push_back(1); else elesFlagMedium->push_back(0); 
      if( isLooseID ) elesFlagLoose->push_back(1); else elesFlagLoose->push_back(0);
      if( isTightID ) elesFlagTight->push_back(1); else elesFlagTight->push_back(0);

      elespfActivity->push_back(pfActivity);
    }
    // add eles to clean from jets
    if( isVetoID && miniIso < maxEleMiniIso_ && ele->pt() > minElePtForElectron2Clean_ ) ele2Clean->push_back(*ele);

      Electron_vetoID->push_back(iPassVetoIDOnly_);//passveto);
      Electron_looseID->push_back(iPassLooseIDOnly_);//passloose);
      Electron_mediumID->push_back(iPassMediumIDOnly_);//passmedium);
      Electron_tightID->push_back(iPassTightIDOnly_);//passtight);

  }


  // determine result before losing ownership of the pointer
  bool result = (doEleVeto_ ? (prod->size() == 0) : true);

  std::unique_ptr<int> nElectrons (new int);

  *nElectrons = prod->size();

  // store in the event
  iEvent.put(std::move(prod));
  iEvent.put(std::move(ele2Clean), "ele2Clean");
  iEvent.put(std::move(elesFlagVeto), "elesFlagVeto");
  iEvent.put(std::move(elesFlagMedium), "elesFlagMedium");
  iEvent.put(std::move(elesFlagTight), "elesFlagTight");
  iEvent.put(std::move(elesFlagLoose), "elesFlagLoose");
  iEvent.put(std::move(elesLVec), "elesLVec");
  iEvent.put(std::move(elesCharge), "elesCharge");
  iEvent.put(std::move(elesMtw), "elesMtw");
  iEvent.put(std::move(elesRelIso), "elesRelIso");
  iEvent.put(std::move(elesisEB), "elesisEB");
  iEvent.put(std::move(elesMiniIso), "elesMiniIso");
  iEvent.put(std::move(elespfActivity), "elespfActivity");
  iEvent.put(std::move(nElectrons), "nElectrons");

  iEvent.put(std::move(Electron_vetoID), "vetoElectronID");
  iEvent.put(std::move(Electron_looseID), "looseElectronID");
  iEvent.put(std::move(Electron_mediumID), "mediumElectronID");
  iEvent.put(std::move(Electron_tightID), "tightElectronID");

  return result;
}

bool prodElectrons::passElectronID(const pat::Electron & ele, const edm::Handle< std::vector<reco::Vertex> > & vertices, const elesIDLevel level) 
{
  // electron ID cuts, updated for Spring15 25ns MC and Run2015C-D data 
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2

  // barrel electrons
  double eb_ieta_cut[4] = {0.0115, 0.011, 0.00998, 0.00998};
  double eb_deta_cut[4] = {0.00749, 0.00477, 0.00311, 0.00308};
  double eb_dphi_cut[4] = {0.228, 0.222, 0.103, 0.0816};
  double eb_hovere_cut[4] = {0.356, 0.298, 0.253, 0.0414};
  double eb_ooeminusoop_cut[4] = {0.299, 0.241, 0.134, 0.0129};
  int eb_misshits_cut[4] = {2, 1, 1, 1};

  // endcap electrons
  double ee_ieta_cut[4] = {0.037, 0.0314, 0.0298, 0.0292};
  double ee_deta_cut[4] = {0.00895, 0.00868, 0.00609, 0.00605};
  double ee_dphi_cut[4] = {0.213, 0.213, 0.045, 0.0394};
  double ee_hovere_cut[4] = {0.211, 0.101, 0.0878, 0.0641};
  double ee_ooeminusoop_cut[4] = {0.15, 0.14, 0.13, 0.0129};
  int ee_misshits_cut[4] = {3, 1, 1, 1};

  // common
  bool reqConvVeto[4] = {true, true, true, true};

  // id variables
  double sigmaIEtaIEta = ele.full5x5_sigmaIetaIeta();
  double dEtaIn        = ele.deltaEtaSuperClusterTrackAtVtx();
  double dPhiIn        = ele.deltaPhiSuperClusterTrackAtVtx();
  double hoe           = ele.hadronicOverEm();
  double ooemoop       = 1e30;
  
  if( ele.ecalEnergy() !=0 && std::isfinite(ele.ecalEnergy()) )
  {
    ooemoop = std::abs(1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy());
  }

  //Impact parameter cuts: differently from previous cut-based electron ID working points, this time the d0 and dz cuts are NOT part of the tuned ID. The reasons for this decision include: efficiency strongly dependent on the physics of the event, measurements wanting to have impact parameter handling different from per-electron d0/dz, relatively low discriminating power of the variables. Instead an average analysis is recommended to use safe highly efficient (when PV is properly found) baseline cuts given in the table below. The d0 and dz are NOT applied in the VID framework, and are left for regular users to cut on explicitly if desired.
  // impact parameter variables
  double eb_d0_cut[4] = {0.05, 0.05, 0.05, 0.05};
  double eb_dz_cut[4] = {0.1, 0.1, 0.1, 0.1};
  double ee_d0_cut[4] = {0.1, 0.1, 0.1, 0.1};
  double ee_dz_cut[4] = {0.2, 0.2, 0.2, 0.2};
  double d0vtx = 0.0;
  double dzvtx = 0.0;
  if (vertices->size() > 0)
  {
    reco::VertexRef vtx(vertices, 0);    
    d0vtx = ele.gsfTrack()->dxy(vtx->position());
    dzvtx = ele.gsfTrack()->dz(vtx->position());
  } 
  else 
  {
    d0vtx = ele.gsfTrack()->dxy();
    dzvtx = ele.gsfTrack()->dz();
  }
  bool dod0dz = dod0dz_, passd0dz_eb = true, passd0dz_ee = true;
  dod0dz ? passd0dz_eb = (eb_d0_cut[level] > fabs(d0vtx)) && (eb_dz_cut[level] > fabs(dzvtx)) : passd0dz_eb = true;
  dod0dz ? passd0dz_ee = (ee_d0_cut[level] > fabs(d0vtx)) && (ee_dz_cut[level] > fabs(dzvtx)) : passd0dz_ee = true;

  // conversion rejection variables
  bool convVeto = ele.passConversionVeto();
  double mHits = ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  
  if (ele.isEB()) 
  {
    return 
         eb_ieta_cut[level] > sigmaIEtaIEta
      && eb_deta_cut[level] > fabs(dEtaIn)
      && eb_dphi_cut[level] > fabs(dPhiIn)
      && eb_hovere_cut[level] > hoe
      && eb_ooeminusoop_cut[level] > fabs(ooemoop)
      && passd0dz_eb
      //&& eb_d0_cut[level] > fabs(d0vtx)
      //&& eb_dz_cut[level] > fabs(dzvtx)
      && (reqConvVeto[level] == convVeto)
      && (eb_misshits_cut[level] >= mHits);
  } 
  else if (ele.isEE()) 
  {
    return 
         ee_ieta_cut[level] > sigmaIEtaIEta
      && ee_deta_cut[level] > fabs(dEtaIn)
      && ee_dphi_cut[level] > fabs(dPhiIn)
      && ee_hovere_cut[level] > hoe
      && ee_ooeminusoop_cut[level] > fabs(ooemoop)
      && passd0dz_ee
      //&& ee_d0_cut[level] > fabs(d0vtx)
      //&& ee_dz_cut[level] > fabs(dzvtx)
      && (reqConvVeto[level] == convVeto)
      && (ee_misshits_cut[level] >= mHits);
  } else return false;

}

bool prodElectrons::passElectronISO(const pat::Electron & ele, const double relIso, const elesIDLevel level)
{
  double eb_relIsoWithEA_cut[4] = {0.175, 0.0994, 0.0695, 0.0588};
  double ee_relIsoWithEA_cut[4] = {0.159, 0.1070, 0.0821, 0.0571}; 
   
  if( ele.isEB() )
  {
    return relIso < eb_relIsoWithEA_cut[level];
  }
  else if (ele.isEE() )
  {
    return relIso < ee_relIsoWithEA_cut[level];
  } 
  else return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(prodElectrons);
