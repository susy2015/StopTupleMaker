
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

#include "DataFormats/METReco/interface/MET.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "TLorentzVector.h"

#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/BoostedDoubleSVTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/IPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "RecoBTag/SecondaryVertex/interface/V0Filter.h"
#include "RecoBTag/ImpactParameter/plugins/IPProducer.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

class prodJets : public edm::EDFilter 
{
 public:

  explicit prodJets(const edm::ParameterSet & iConfig);
  ~prodJets();

 private:

  virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);
  void compute(const reco::Jet * jet, bool isReco, double& totalMult_, double& ptD_, double& axis1_, double& axis2_);

  edm::InputTag jetSrc_, jetOtherSrc_;
  // All have to be pat::Jet, otherwise cannot get b-tagging information!
  edm::Handle<std::vector<pat::Jet> > jets, otherjets; 
  std::string bTagKeyString_;
  edm::InputTag vtxSrc_;
  edm::InputTag metSrc_;
  bool isPatJet;
  bool debug_;

  bool isData_;

  double jetPtCut_miniAOD_, genMatch_dR_;
  double relPt_for_xCheck_, dR_for_xCheck_;

  edm::EDGetTokenT<std::vector<pat::Jet> >JetTok_;
  edm::EDGetTokenT<std::vector<pat::Jet> >OtherJetsTok_;
  edm::EDGetTokenT<std::vector<int> > W_EmuVec_Tok_;
  edm::EDGetTokenT<std::vector<int> >W_TauVec_Tok_;
  edm::EDGetTokenT<std::vector<int> >W_Tau_EmuVec_Tok_;
  edm::EDGetTokenT<std::vector<int> >W_Tau_ProngsVec_Tok_;
  edm::EDGetTokenT<std::vector<int> >W_Tau_NuVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >GenDecayLVec_Tok_;
  edm::EDGetTokenT<std::vector<int> >GenDecayMomRefVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >EleLVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >MuLVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >TrksForIsoVetolVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >LooseIsoTrksVec_Tok_;
  edm::EDGetTokenT< std::vector<reco::Vertex> >VtxTok_;
  edm::EDGetTokenT<std::vector<pat::Jet>> PuppiJetsSrc_Tok_;
  edm::EDGetTokenT<std::vector<pat::Jet>> PuppiSubJetsSrc_Tok_;
  edm::EDGetTokenT<std::vector<pat::Jet>> Ak8Jets_Tok_;
  edm::EDGetTokenT<std::vector<pat::Jet>> Ak8SubJets_Tok_;

  edm::InputTag W_emuVec_Src_, W_tauVec_Src_, W_tau_emuVec_Src_, W_tau_prongsVec_Src_, W_tau_nuVec_Src_;
  edm::Handle<std::vector<int> > W_emuVec_, W_tauVec_, W_tau_emuVec_, W_tau_prongsVec_, W_tau_nuVec_;

  edm::InputTag genDecayLVec_Src_;
  edm::Handle<std::vector<TLorentzVector> > genDecayLVec_;

  edm::InputTag genDecayMomRefVec_Src_;
  edm::Handle<std::vector<int> > genDecayMomRefVec_;

  edm::InputTag eleLVec_Src_, muLVec_Src_;
  edm::Handle<std::vector<TLorentzVector> > eleLVec_, muLVec_;

  edm::InputTag trksForIsoVetoLVec_Src_, looseisoTrksLVec_Src_;
  edm::Handle<std::vector<TLorentzVector> > trksForIsoVetoLVec_, looseisoTrksLVec_;
  double deltaRcon_;

  //PUPPI sources
  edm::InputTag puppiJetsSrc_, puppiSubJetsSrc_;
  edm::Handle<std::vector<pat::Jet> > puppiJets;
  edm::Handle<std::vector<pat::Jet> > puppiSubJets; 
 
  std::string jetType_;
  std::string qgTaggerKey_;
  std::string NjettinessAK8Puppi_label_;
  std::string ak8PFJetsPuppi_label_;

  std::string jetPBJetTags_;
  std::string jetPNegBJetTags_;
  std::string jetPPosBJetTags_;

  std::string jetBPBJetTags_;
  std::string jetBPNegBJetTags_;
  std::string jetBPPosBJetTags_;

  std::string deepCSVBJetTags_;
  std::string deepCSVNegBJetTags_;
  std::string deepCSVPosBJetTags_;

  std::string deepFlavorBJetTags_;

 std::string combinedSVBJetTags_;
  std::string combinedSVNegBJetTags_;
  std::string combinedSVPosBJetTags_;

  std::string combinedIVFSVBJetTags_;
  std::string combinedIVFSVPosBJetTags_;
  std::string combinedIVFSVNegBJetTags_;

  std::string simpleSVHighEffBJetTags_;
  std::string simpleSVNegHighEffBJetTags_;
  std::string simpleSVHighPurBJetTags_;
  std::string simpleSVNegHighPurBJetTags_;

  std::string softPFMuonBJetTags_;
  std::string softPFMuonNegBJetTags_;
  std::string softPFMuonPosBJetTags_;

  std::string softPFElectronBJetTags_;
  std::string softPFElectronNegBJetTags_;
  std::string softPFElectronPosBJetTags_;

  std::string doubleSVBJetTags_;

  std::string cMVABJetTags_;
  std::string cMVAv2BJetTags_;
  std::string cMVAv2NegBJetTags_;
  std::string cMVAv2PosBJetTags_;

  std::string   CvsBCJetTags_;
  std::string   CvsBNegCJetTags_;
  std::string   CvsBPosCJetTags_;
  std::string   CvsLCJetTags_;
  std::string   CvsLNegCJetTags_;
  std::string   CvsLPosCJetTags_;
};

void prodJets::compute(const reco::Jet * jet, bool isReco, double& totalMult_, double& ptD_, double& axis1_, double& axis2_)
{
    totalMult_ = 0;
    ptD_       = 0;
    axis1_     = 0;
    axis2_     = 0;
    if(jet->numberOfDaughters() == 0) return;
    float sum_weight    = 0.0;
    float sum_dEta      = 0.0;
    float sum_dPhi      = 0.0;
    float sum_dEta2     = 0.0;
    float sum_dPhi2     = 0.0;
    float sum_dEta_dPhi = 0.0;
    float sum_pt        = 0.0;
    bool useQC          = false; // useQualityCuts; hard-coded for now to mimic what jetMet does in 731

    // loop over the jet constituents
    // (packed candidate situation)
    for(auto part : jet->getJetConstituentsQuick()) {
        if(part->charge()){ // charged particles
          if(isReco) {
               auto p = dynamic_cast<const pat::PackedCandidate*>(part);
               if(!p){
                      try { throw; }
                      catch(...) {
                          std::cout << "ERROR: QGTagging variables cannot be computed for these jets!" << std::endl
                              << "       See QuauarGluonTaggingVaiables::compute()"              << std::endl;
                      } // catch(...)
               } // !p
               if(!( p->fromPV() > 1 && p->trackHighPurity() )) continue;
               if(useQC) {
                   // currently hard-coded to false above
                   // this isn't stored for packedCandidates, so will need fix if useQC is changed to true
                   if( p->dzError()==0 || p->dxyError()==0 ) continue;
                   if( (p->dz()*p->dz() )  / (p->dzError()*p->dzError() ) > 25. ) continue;
                    if( (p->dxy()*p->dxy()) / (p->dxyError()*p->dxyError()) < 25. ) ++totalMult_; // this cut only applies to multiplicity
                                    } else ++totalMult_;
          } else ++totalMult_;
        } else { // neutral particles
            if(part->pt() < 1.0) continue;
            ++totalMult_;
        } // charged, neutral particles
        float dEta   = part->eta() - jet->eta();
        float dPhi   = reco::deltaPhi(part->phi(), jet->phi());
        float partPt = part->pt();
        float weight = partPt*partPt;

        sum_weight    += weight;
        sum_pt        += partPt;
        sum_dEta      += dEta      * weight;
        sum_dPhi      += dPhi      * weight;
        sum_dEta2     += dEta*dEta * weight;
        sum_dEta_dPhi += dEta*dPhi * weight;
        sum_dPhi2     += dPhi*dPhi * weight;
    } // jet->getJetConstituentsQuick()
// calculate axis2 and ptD
    float a = 0.0;
    float b = 0.0;
    float c = 0.0;
    float ave_dEta  = 0.0;
    float ave_dPhi  = 0.0;
    float ave_dEta2 = 0.0;
    float ave_dPhi2 = 0.0;

    if(sum_weight > 0){
    ptD_ = sqrt(sum_weight)/sum_pt;
    ave_dEta  = sum_dEta  / sum_weight;
    ave_dPhi  = sum_dPhi  / sum_weight;
    ave_dEta2 = sum_dEta2 / sum_weight;
    ave_dPhi2 = sum_dPhi2 / sum_weight;
    a = ave_dEta2 - ave_dEta*ave_dEta;
    b = ave_dPhi2 - ave_dPhi*ave_dPhi;
    c = -(sum_dEta_dPhi/sum_weight - ave_dEta*ave_dPhi);
    } else ptD_ = 0;    
    float delta = sqrt(fabs( (a-b)*(a-b) + 4*c*c ));
    if(a+b-delta > 0) axis2_ = sqrt(0.5*(a+b-delta));
        else              axis2_ = 0.0;
           if(a+b+delta > 0) axis1_ = sqrt(0.5*(a+b+delta));
            else              axis1_ = 0.0;
}

prodJets::prodJets(const edm::ParameterSet & iConfig) 
{
  isData_ = true;

  jetSrc_      = iConfig.getParameter<edm::InputTag>("jetSrc");
  jetOtherSrc_ = iConfig.getParameter<edm::InputTag>("jetOtherSrc");
  vtxSrc_      = iConfig.getParameter<edm::InputTag>("vtxSrc");
  //metSrc_      = iConfig.getParameter<edm::InputTag>("metSrc");
  bTagKeyString_ = iConfig.getParameter<std::string>("bTagKeyString");

  jetPBJetTags_        = iConfig.getParameter<std::string>("jetPBJetTags");
  jetPNegBJetTags_     = iConfig.getParameter<std::string>("jetPNegBJetTags");
  jetPPosBJetTags_     = iConfig.getParameter<std::string>("jetPPosBJetTags");

  jetBPBJetTags_        = iConfig.getParameter<std::string>("jetBPBJetTags");
  jetBPNegBJetTags_     = iConfig.getParameter<std::string>("jetBPNegBJetTags");
  jetBPPosBJetTags_     = iConfig.getParameter<std::string>("jetBPPosBJetTags");

  deepCSVBJetTags_    = iConfig.getParameter<std::string>("deepCSVBJetTags");
  deepCSVNegBJetTags_ = iConfig.getParameter<std::string>("deepCSVNegBJetTags");
  deepCSVPosBJetTags_ = iConfig.getParameter<std::string>("deepCSVPosBJetTags");

  deepFlavorBJetTags_    = iConfig.getParameter<std::string>("deepFlavorBJetTags");

  combinedSVBJetTags_     = iConfig.getParameter<std::string>("combinedSVBJetTags");
  combinedSVNegBJetTags_  = iConfig.getParameter<std::string>("combinedSVNegBJetTags");
  combinedSVPosBJetTags_  = iConfig.getParameter<std::string>("combinedSVPosBJetTags");

  combinedIVFSVBJetTags_      = iConfig.getParameter<std::string>("combinedIVFSVBJetTags");
  combinedIVFSVPosBJetTags_   = iConfig.getParameter<std::string>("combinedIVFSVPosBJetTags");
  combinedIVFSVNegBJetTags_   = iConfig.getParameter<std::string>("combinedIVFSVNegBJetTags");

  simpleSVHighEffBJetTags_      = iConfig.getParameter<std::string>("simpleSVHighEffBJetTags");
  simpleSVNegHighEffBJetTags_   = iConfig.getParameter<std::string>("simpleSVNegHighEffBJetTags");
  simpleSVHighPurBJetTags_      = iConfig.getParameter<std::string>("simpleSVHighPurBJetTags");
  simpleSVNegHighPurBJetTags_   = iConfig.getParameter<std::string>("simpleSVNegHighPurBJetTags");

  combinedIVFSVBJetTags_      = iConfig.getParameter<std::string>("combinedIVFSVBJetTags");
  combinedIVFSVPosBJetTags_   = iConfig.getParameter<std::string>("combinedIVFSVPosBJetTags");
  combinedIVFSVNegBJetTags_   = iConfig.getParameter<std::string>("combinedIVFSVNegBJetTags");

  softPFMuonBJetTags_       = iConfig.getParameter<std::string>("softPFMuonBJetTags");
  softPFMuonNegBJetTags_    = iConfig.getParameter<std::string>("softPFMuonNegBJetTags");
  softPFMuonPosBJetTags_    = iConfig.getParameter<std::string>("softPFMuonPosBJetTags");

  softPFElectronBJetTags_       = iConfig.getParameter<std::string>("softPFElectronBJetTags");
  softPFElectronNegBJetTags_    = iConfig.getParameter<std::string>("softPFElectronNegBJetTags");
  softPFElectronPosBJetTags_    = iConfig.getParameter<std::string>("softPFElectronPosBJetTags");

  doubleSVBJetTags_ = iConfig.getParameter<std::string>("doubleSVBJetTags");

  cMVABJetTags_ = iConfig.getParameter<std::string>("cMVABJetTags");
  cMVAv2BJetTags_ = iConfig.getParameter<std::string>("cMVAv2BJetTags");
  cMVAv2NegBJetTags_ = iConfig.getParameter<std::string>("cMVAv2NegBJetTags");
  cMVAv2PosBJetTags_ = iConfig.getParameter<std::string>("cMVAv2PosBJetTags");

  CvsBCJetTags_             = iConfig.getParameter<std::string>("CvsBCJetTags");
  CvsBNegCJetTags_             = iConfig.getParameter<std::string>("CvsBNegCJetTags");
  CvsBPosCJetTags_             = iConfig.getParameter<std::string>("CvsBPosCJetTags");
  CvsLCJetTags_             = iConfig.getParameter<std::string>("CvsLCJetTags");
  CvsLNegCJetTags_             = iConfig.getParameter<std::string>("CvsLNegCJetTags");
  CvsLPosCJetTags_             = iConfig.getParameter<std::string>("CvsLPosCJetTags");

  deepFlavorBJetTags_    = iConfig.getParameter<std::string>("deepFlavorBJetTags");

  debug_       = iConfig.getParameter<bool>("debug");

  jetPtCut_miniAOD_ = iConfig.getUntrackedParameter<double>("jetPtCut_miniAOD", 10);
  genMatch_dR_ = iConfig.getUntrackedParameter<double>("genMatch_dR", 1.0);
  dR_for_xCheck_ = iConfig.getUntrackedParameter<double>("dR_for_xCheck", 0.2);
  relPt_for_xCheck_ = iConfig.getUntrackedParameter<double>("relPt_for_xCheck", 1e-2);

  W_emuVec_Src_ = iConfig.getParameter<edm::InputTag>("W_emuVec");
  W_tauVec_Src_ = iConfig.getParameter<edm::InputTag>("W_tauVec");
  W_tau_emuVec_Src_ = iConfig.getParameter<edm::InputTag>("W_tau_emuVec");
  W_tau_prongsVec_Src_ = iConfig.getParameter<edm::InputTag>("W_tau_prongsVec");
  W_tau_nuVec_Src_ = iConfig.getParameter<edm::InputTag>("W_tau_nuVec");

  genDecayLVec_Src_ = iConfig.getParameter<edm::InputTag>("genDecayLVec");

  genDecayMomRefVec_Src_ = iConfig.getParameter<edm::InputTag>("genDecayMomRefVec");

  eleLVec_Src_ = iConfig.getParameter<edm::InputTag>("eleLVec");
  muLVec_Src_ = iConfig.getParameter<edm::InputTag>("muLVec");
  
  trksForIsoVetoLVec_Src_ = iConfig.getParameter<edm::InputTag>("trksForIsoVetoLVec");
  looseisoTrksLVec_Src_ = iConfig.getParameter<edm::InputTag>("looseisoTrksLVec");

  deltaRcon_ = iConfig.getUntrackedParameter<double>("deltaRcon", 0.01);

  jetType_ = iConfig.getParameter<std::string>("jetType");

  qgTaggerKey_ = iConfig.getParameter<std::string>("qgTaggerKey");
  
  //Puppi source
  puppiJetsSrc_ = iConfig.getParameter<edm::InputTag>("puppiJetsSrc");
  puppiSubJetsSrc_ = iConfig.getParameter<edm::InputTag>("puppiSubJetsSrc");

  NjettinessAK8Puppi_label_ = iConfig.getParameter<std::string>("NjettinessAK8Puppi_label");
  ak8PFJetsPuppi_label_ = iConfig.getParameter<std::string>("ak8PFJetsPuppi_label");


  JetTok_ = consumes<std::vector<pat::Jet> >(jetSrc_);
  OtherJetsTok_ = consumes<std::vector<pat::Jet> >(jetOtherSrc_);
  W_EmuVec_Tok_=consumes<std::vector<int> >(W_emuVec_Src_);
  W_TauVec_Tok_=consumes<std::vector<int> >(W_tauVec_Src_);
  W_Tau_EmuVec_Tok_=consumes<std::vector<int> >(W_tau_emuVec_Src_);
  W_Tau_ProngsVec_Tok_ = consumes<std::vector<int> >(W_tau_prongsVec_Src_);
  W_Tau_NuVec_Tok_ = consumes<std::vector<int> >(W_tau_nuVec_Src_);
  GenDecayLVec_Tok_=consumes<std::vector<TLorentzVector> >(genDecayLVec_Src_);
  GenDecayMomRefVec_Tok_=consumes<std::vector<int> >(genDecayMomRefVec_Src_);
  EleLVec_Tok_=consumes<std::vector<TLorentzVector> >(eleLVec_Src_);
  MuLVec_Tok_=consumes<std::vector<TLorentzVector> >(muLVec_Src_);
  TrksForIsoVetolVec_Tok_=consumes<std::vector<TLorentzVector> >(trksForIsoVetoLVec_Src_);
  LooseIsoTrksVec_Tok_=consumes<std::vector<TLorentzVector> >(looseisoTrksLVec_Src_);
  VtxTok_=consumes< std::vector<reco::Vertex> >(vtxSrc_);
  PuppiJetsSrc_Tok_ = consumes<std::vector<pat::Jet>>(puppiJetsSrc_);
  PuppiSubJetsSrc_Tok_ = consumes<std::vector<pat::Jet>>(puppiSubJetsSrc_);

  //produces<std::vector<pat::Jet> >("");
  produces<std::vector<TLorentzVector> >("jetsLVec");
  produces<std::vector<int> >("recoJetsFlavor");
  produces<std::vector<float> >("recoJetsBtag");
  produces<std::vector<float> >("recoJetsCharge");
  produces<std::vector<float> >("recoJetsJecUnc");
  produces<std::vector<float> >("recoJetsJecScaleRawToFull");
  produces<int>("nJets");
  produces<std::vector<float> >("qgLikelihood");
  produces<std::vector<float> >("qgPtD");
  produces<std::vector<float> >("qgAxis2");
  produces<std::vector<float> >("qgAxis1");
  produces<std::vector<int> >("qgMult");

  //produce variables needed for Lost Lepton study, added by hua.wei@cern.ch
  produces<std::vector<float> >("recoJetschargedHadronEnergyFraction");
  produces<std::vector<float> >("recoJetschargedEmEnergyFraction");
  produces<std::vector<float> >("recoJetsneutralEmEnergyFraction");

  produces<std::vector<float> >("recoJetsmuonEnergyFraction");

  produces<std::vector<int> >("muMatchedJetIdx");
  produces<std::vector<int> >("eleMatchedJetIdx");

  produces<std::vector<int> >("trksForIsoVetoMatchedJetIdx");
  produces<std::vector<int> >("looseisoTrksMatchedJetIdx");

  //PUPPI
  produces<std::vector<TLorentzVector> >("puppiJetsLVec");
  produces<std::vector<TLorentzVector> >("puppiSubJetsLVec");
  produces<std::vector<float> >("puppisoftDropMass");
  produces<std::vector<float> >("puppitau1");
  produces<std::vector<float> >("puppitau2");
  produces<std::vector<float> >("puppitau3");
  produces<std::vector<float> >("puppiSubJetsBdisc");

produces<std::vector<float> >("JetProba");
  produces<std::vector<float> >("JetProbaN");
  produces<std::vector<float> >("JetProbaP");
  produces<std::vector<float> >("JetBprob");
  produces<std::vector<float> >("JetBprobP");
  produces<std::vector<float> >("JetBprobN");

  produces<std::vector<float> >("CombinedSvtx");
  produces<std::vector<float> >("CombinedSvtxN");
  produces<std::vector<float> >("CombinedSvtxP");

produces<std::vector<float> >("DeepCSVb");
  produces<std::vector<float> >("DeepCSVc");
  produces<std::vector<float> >("DeepCSVl");
  produces<std::vector<float> >("DeepCSVbb");
  produces<std::vector<float> >("DeepCSVcc");

  produces<std::vector<float> >("DeepCSVbN");
  produces<std::vector<float> >("DeepCSVcN");
  produces<std::vector<float> >("DeepCSVlN");
  produces<std::vector<float> >("DeepCSVbbN");
  produces<std::vector<float> >("DeepCSVccN");

  produces<std::vector<float> >("DeepCSVbP");
  produces<std::vector<float> >("DeepCSVcP");
  produces<std::vector<float> >("DeepCSVlP");
  produces<std::vector<float> >("DeepCSVbbP");
  produces<std::vector<float> >("DeepCSVccP");

  produces<std::vector<float> >("DeepFlavorb");
  produces<std::vector<float> >("DeepFlavorbb");
  produces<std::vector<float> >("DeepFlavorlepb");
  produces<std::vector<float> >("DeepFlavorc");
  produces<std::vector<float> >("DeepFlavoruds");
  produces<std::vector<float> >("DeepFlavorg");

  produces<std::vector<float> >("CombinedIVF");
  produces<std::vector<float> >("CombinedIVFP");
  produces<std::vector<float> >("CombinedIVFN");

  produces<std::vector<float> >("Svtx");
  produces<std::vector<float> >("SvtxN");
  produces<std::vector<float> >("SvtxHP");
  produces<std::vector<float> >("SvtxNHP");

  produces<std::vector<float> >("SoftM");
  produces<std::vector<float> >("SoftMN");
  produces<std::vector<float> >("SoftMP");

  produces<std::vector<float> >("SoftE");
  produces<std::vector<float> >("SoftEN");
  produces<std::vector<float> >("SoftEP");

  produces<std::vector<float> >("DoubleSV");

  produces<std::vector<float> >("cMVA");
  produces<std::vector<float> >("cMVAv2");
  produces<std::vector<float> >("cMVAv2Neg");
  produces<std::vector<float> >("cMVAv2Pos");

  produces<std::vector<float> >("CvsB");
  produces<std::vector<float> >("CvsBNeg");
  produces<std::vector<float> >("CvsBPos");
  produces<std::vector<float> >("CvsL");
  produces<std::vector<float> >("CvsLNeg");
  produces<std::vector<float> >("CvsLPos");

}


prodJets::~prodJets() 
{
}


bool prodJets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  if( !iEvent.isRealData() ) isData_ = false;

  iEvent.getByToken(JetTok_, jets);

  //get the JEC uncertainties
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jetType_, JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  std::unique_ptr<JetCorrectionUncertainty> jecUnc( new JetCorrectionUncertainty(JetCorPar) );

  if( !isData_ )
  {
    iEvent.getByToken(OtherJetsTok_, otherjets);
    iEvent.getByToken(W_EmuVec_Tok_, W_emuVec_);
    iEvent.getByToken(W_TauVec_Tok_, W_tauVec_);
    iEvent.getByToken(W_Tau_EmuVec_Tok_, W_tau_emuVec_);
    iEvent.getByToken(W_Tau_ProngsVec_Tok_, W_tau_prongsVec_);
    iEvent.getByToken(W_Tau_NuVec_Tok_, W_tau_nuVec_);

    iEvent.getByToken(GenDecayLVec_Tok_, genDecayLVec_);
    iEvent.getByToken(GenDecayMomRefVec_Tok_, genDecayMomRefVec_);
  }

  iEvent.getByToken(EleLVec_Tok_, eleLVec_);
  iEvent.getByToken(MuLVec_Tok_, muLVec_);

  iEvent.getByToken(TrksForIsoVetolVec_Tok_, trksForIsoVetoLVec_);
  iEvent.getByToken(LooseIsoTrksVec_Tok_,looseisoTrksLVec_);

  //read in the objects
  edm::Handle< std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(VtxTok_, vertices);
  //reco::Vertex::Point vtxpos = (vertices->size() > 0 ? (*vertices)[0].position() : reco::Vertex::Point());
  //edm::Handle<edm::View<reco::MET> > met;
  //iEvent.getByLabel(metSrc_, met);

  std::vector<pat::Jet> extJets = (*jets);
  //PUPPI
  iEvent.getByToken(PuppiJetsSrc_Tok_, puppiJets);
  iEvent.getByToken(PuppiSubJetsSrc_Tok_, puppiSubJets);

  //check which ones to keep
  //std::unique_ptr<std::vector<pat::Jet> > prod(new std::vector<pat::Jet>());
  std::unique_ptr<std::vector<TLorentzVector> > jetsLVec(new std::vector<TLorentzVector>());
  std::unique_ptr<std::vector<int> > recoJetsFlavor(new std::vector<int>());
  std::unique_ptr<std::vector<float> > recoJetsBtag(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsCharge(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsJecUnc(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsJecScaleRawToFull(new std::vector<float>());

 std::unique_ptr<std::vector<float> > JetProba(new std::vector<float>());
  std::unique_ptr<std::vector<float> > JetProbaN(new std::vector<float>());
  std::unique_ptr<std::vector<float> > JetProbaP(new std::vector<float>());
  std::unique_ptr<std::vector<float> > JetBprob(new std::vector<float>());
  std::unique_ptr<std::vector<float> > JetBprobP(new std::vector<float>());
  std::unique_ptr<std::vector<float> > JetBprobN(new std::vector<float>());

  std::unique_ptr<std::vector<float> >CombinedSvtx(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CombinedSvtxN(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CombinedSvtxP(new std::vector<float>());

  std::unique_ptr<std::vector<float> > DeepCSVb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVc(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVl(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVbb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVcc(new std::vector<float>());

  std::unique_ptr<std::vector<float> > DeepCSVbN(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVcN(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVlN(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVbbN(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVccN(new std::vector<float>());

  std::unique_ptr<std::vector<float> > DeepCSVbP(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVcP(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVlP(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVbbP(new std::vector<float>());
std::unique_ptr<std::vector<float> > DeepCSVccP(new std::vector<float>());

  std::unique_ptr<std::vector<float> > DeepFlavorb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorbb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorlepb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorc(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavoruds(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorg(new std::vector<float>());

  std::unique_ptr<std::vector<float> >CombinedIVF(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CombinedIVFN(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CombinedIVFP(new std::vector<float>());

  std::unique_ptr<std::vector<float> >Svtx(new std::vector<float>());
  std::unique_ptr<std::vector<float> >SvtxN(new std::vector<float>());
  std::unique_ptr<std::vector<float> >SvtxHP(new std::vector<float>());
  std::unique_ptr<std::vector<float> >SvtxNHP(new std::vector<float>());

  std::unique_ptr<std::vector<float> >SoftM(new std::vector<float>());
  std::unique_ptr<std::vector<float> >SoftMN(new std::vector<float>());
  std::unique_ptr<std::vector<float> >SoftMP(new std::vector<float>());

  std::unique_ptr<std::vector<float> >SoftE(new std::vector<float>());
  std::unique_ptr<std::vector<float> >SoftEN(new std::vector<float>());
  std::unique_ptr<std::vector<float> >SoftEP(new std::vector<float>());

  std::unique_ptr<std::vector<float> >DoubleSV(new std::vector<float>());

  std::unique_ptr<std::vector<float> >cMVA(new std::vector<float>());
  std::unique_ptr<std::vector<float> >cMVAv2(new std::vector<float>());
  std::unique_ptr<std::vector<float> >cMVAv2Neg(new std::vector<float>());
  std::unique_ptr<std::vector<float> >cMVAv2Pos(new std::vector<float>());

  std::unique_ptr<std::vector<float> >CvsB(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CvsBNeg(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CvsBPos(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CvsL(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CvsLNeg(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CvsLPos(new std::vector<float>());

  std::unique_ptr<std::vector<float> > qgLikelihood(new std::vector<float>());
  std::unique_ptr<std::vector<float> > qgPtD(new std::vector<float>());
  std::unique_ptr<std::vector<float> > qgAxis2(new std::vector<float>());
  std::unique_ptr<std::vector<float> > qgAxis1(new std::vector<float>());
  std::unique_ptr<std::vector<int> > qgMult(new std::vector<int>());

  std::unique_ptr<std::vector<float> > recoJetschargedHadronEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetschargedEmEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsneutralEmEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsmuonEnergyFraction(new std::vector<float>());

  std::unique_ptr<std::vector<int> > muMatchedJetIdx(new std::vector<int>(muLVec_->size(), -1));
  std::unique_ptr<std::vector<int> > eleMatchedJetIdx(new std::vector<int>(eleLVec_->size(), -1));

  std::unique_ptr<std::vector<int> > trksForIsoVetoMatchedJetIdx(new std::vector<int>(trksForIsoVetoLVec_->size(), -1));
  std::unique_ptr<std::vector<int> > looseisoTrksMatchedJetIdx(new std::vector<int>(looseisoTrksLVec_->size(), -1));

  //PUPPI
  std::unique_ptr<std::vector<TLorentzVector> > puppiJetsLVec(new std::vector<TLorentzVector>());
  std::unique_ptr<std::vector<TLorentzVector> > puppiSubJetsLVec(new std::vector<TLorentzVector>());
  std::unique_ptr<std::vector<float> > puppiSubJetsBdisc(new std::vector<float>());
  std::unique_ptr<std::vector<float> > puppisoftDropMass(new std::vector<float>());
  std::unique_ptr<std::vector<float> > puppitau1(new std::vector<float>());
  std::unique_ptr<std::vector<float> > puppitau2(new std::vector<float>());
  std::unique_ptr<std::vector<float> > puppitau3(new std::vector<float>());



  if( !isData_ ){
     int cntJetPassPtCut = 0;
     for(unsigned int io=0; io < otherjets->size(); io++){
        const double otjet_pt = otherjets->at(io).pt(), otjet_eta = otherjets->at(io).eta(), otjet_phi = otherjets->at(io).phi();
        TLorentzVector perLVec; perLVec.SetPtEtaPhiE(otjet_pt, otjet_eta, otjet_phi, otherjets->at(io).energy());
        int cntFound = 0, matchedIdx = -1;
        double minDR = 999.0;
        for(unsigned int ij=0; ij< jets->size(); ij++){
           const double jet_eta = jets->at(ij).eta(), jet_phi = jets->at(ij).phi();
           const double dR = reco::deltaR(otjet_eta, otjet_phi, jet_eta, jet_phi);
           if( minDR > dR ){
              minDR = dR; matchedIdx = ij;
           }
        }
      if( matchedIdx != -1 )
      {
        if( minDR < dR_for_xCheck_ && std::abs(otjet_pt - jets->at(matchedIdx).pt())/jets->at(matchedIdx).pt() < relPt_for_xCheck_ ){ cntFound ++; }
      }
      if( otjet_pt >= jetPtCut_miniAOD_ )
      {
        cntJetPassPtCut ++;
        if( cntFound != 1 && debug_ )
        {
          std::cout<<"WARNING ... jet mis-matching between otherjets and jets for pt > "<<jetPtCut_miniAOD_<<"  matchedIdx : "<<matchedIdx<<"  cntFound : "<<cntFound<<"  minDR : "<<minDR<<std::endl;
          std::cout<<"otjet_pt : "<<otjet_pt<<"  otjet_eta : "<<otjet_eta<<"  otjet_phi : "<<otjet_phi<<std::endl;
          if( matchedIdx != -1 ) std::cout<<"  jet_pt : "<<jets->at(matchedIdx).pt()<<"    jet_eta : "<<jets->at(matchedIdx).eta()<<"    jet_phi : "<<jets->at(matchedIdx).phi()<<std::endl;
        }
      }
      else
      {
        if( cntFound && debug_ )
        {
          std::cout<<"WARNING ... otjet with pt : "<<otjet_pt<<"  matching to one of the jets for pt > "<<jetPtCut_miniAOD_<<" ?!"<<std::endl;
          std::cout<<"otjet_pt : "<<otjet_pt<<"  otjet_eta : "<<otjet_eta<<"  otjet_phi : "<<otjet_phi<<std::endl;
          std::cout<<"  jet_pt : "<<jets->at(matchedIdx).pt()<<"    jet_eta : "<<jets->at(matchedIdx).eta()<<"    jet_phi : "<<jets->at(matchedIdx).phi()<<std::endl;
        }
        else
        {
          int cntgenMatch = 0;
          for(unsigned int ig=0; ig<W_emuVec_->size(); ig++)
          {
            int perIdx = W_emuVec_->at(ig);
            TLorentzVector genLVec = genDecayLVec_->at(perIdx);
            const double perdeltaR = perLVec.DeltaR(genLVec);
            if( perdeltaR < genMatch_dR_ ) cntgenMatch ++;
          }
          for(unsigned int ig=0; ig<W_tauVec_->size(); ig++)
          {
            int perIdx = W_tauVec_->at(ig);
            TLorentzVector genLVec = genDecayLVec_->at(perIdx);
            const double perdeltaR = perLVec.DeltaR(genLVec);
            if( perdeltaR < genMatch_dR_ ) cntgenMatch ++;
          }
          for(unsigned int ig=0; ig<W_tau_emuVec_->size(); ig++)
          {
            int perIdx = W_tau_emuVec_->at(ig);
            TLorentzVector genLVec = genDecayLVec_->at(perIdx);
            const double perdeltaR = perLVec.DeltaR(genLVec);
            if( perdeltaR < genMatch_dR_ ) cntgenMatch ++;
          }
          for(unsigned int ig=0; ig<W_tau_prongsVec_->size(); ig++)
          {
            int perIdx = W_tau_prongsVec_->at(ig);
            TLorentzVector genLVec = genDecayLVec_->at(perIdx);
            const double perdeltaR = perLVec.DeltaR(genLVec);
            if( perdeltaR < genMatch_dR_ ) cntgenMatch ++;
          }
          for(unsigned int ig=0; ig<W_tauVec_->size(); ig++)
          {
            int perIdx = W_tauVec_->at(ig);
            TLorentzVector genLVec = genDecayLVec_->at(perIdx);
            for(unsigned int in=0; in<W_tau_nuVec_->size(); in++)
            {
              int perJdx = W_tau_nuVec_->at(in);
              TLorentzVector gennuLVec = genDecayLVec_->at(perJdx);
              int momIdx = perJdx;
              bool isFound = false;
              while( momIdx != -1 )
              {
                momIdx = genDecayMomRefVec_->at(momIdx);
                if( momIdx == perIdx ){ isFound = true; break; }
              }
              if( isFound ) genLVec -= gennuLVec;
            }
            const double perdeltaR = perLVec.DeltaR(genLVec);
            if( perdeltaR < genMatch_dR_ ) cntgenMatch ++;
          }
          for(unsigned int im=0; im<muLVec_->size(); im++)
          {
            const double perdeltaR = perLVec.DeltaR(muLVec_->at(im));
            if( perdeltaR < genMatch_dR_ ) cntgenMatch ++;
          }
          for(unsigned int ie=0; ie<eleLVec_->size(); ie++)
          {
            const double perdeltaR = perLVec.DeltaR(eleLVec_->at(ie));
            if( perdeltaR < genMatch_dR_ ) cntgenMatch ++;
          }
          if( cntgenMatch ){ extJets.push_back(otherjets->at(io)); }
        }
      }
    } 
   
    if( cntJetPassPtCut != (int)jets->size() && debug_ ) std::cout<<"WARNING ... cntJetPassPtCut : "<<cntJetPassPtCut<<"  NOT EQUAL jets->size : "<<jets->size()<<std::endl;
    if( (int)jets->size() >= 4 && std::abs(1.0*cntJetPassPtCut - 1.0*jets->size())/(1.0*jets->size()) > 0.1 && debug_ )
    {
      std::cout<<"\nWARNING ... cntJetPassPtCut : "<<cntJetPassPtCut<<"  slimmedJets.size : "<<jets->size()<<std::endl;
      std::cout<<"Please checking if global tag used for re-clustering the jets is the same as used to produce the miniAOD!"<<std::endl<<std::endl;
    }
  }
  //PUPPI

  for(unsigned int ip = 0; ip < puppiJets->size(); ip++){
     TLorentzVector perPuppiJetLVec;
     perPuppiJetLVec.SetPtEtaPhiE( puppiJets->at(ip).pt(), puppiJets->at(ip).eta(), puppiJets->at(ip).phi(), puppiJets->at(ip).energy() );
     puppiJetsLVec->push_back(perPuppiJetLVec);

     float puppi_tau1_uf         = puppiJets->at(ip).userFloat(NjettinessAK8Puppi_label_+":tau1");
     puppitau1->push_back(puppi_tau1_uf);
     float puppi_tau2_uf         = puppiJets->at(ip).userFloat(NjettinessAK8Puppi_label_+":tau2");
     puppitau2->push_back(puppi_tau2_uf);
     float puppi_tau3_uf         = puppiJets->at(ip).userFloat(NjettinessAK8Puppi_label_+":tau3");
     puppitau3->push_back(puppi_tau3_uf);

     float puppisoftDropMass_uf = puppiJets->at(ip).userFloat(ak8PFJetsPuppi_label_+"SoftDropMass");
     puppisoftDropMass->push_back(puppisoftDropMass_uf);
  }

  for(unsigned int ip = 0; ip < puppiSubJets->size(); ip++){
     // The subjet collection is a collection of softdropped jets, and you have to access the subjets from there
     // Most of the time there are two daughters, sometimes there is only one
     // Not all of the puppiJets have a matching puppiSubJets jet

     TLorentzVector perPuppiSubJetLVec;
     for(unsigned int id=0; id<puppiSubJets->at(ip).numberOfDaughters(); ++id){
       perPuppiSubJetLVec.SetPtEtaPhiE( puppiSubJets->at(ip).daughter(id)->pt(), 
					puppiSubJets->at(ip).daughter(id)->eta(), 
					puppiSubJets->at(ip).daughter(id)->phi(), 
					puppiSubJets->at(ip).daughter(id)->energy() );
       puppiSubJetsLVec->push_back(perPuppiSubJetLVec);
       // btag info is not available for reco::Candidates, so must cast it to a pat::Jet
       puppiSubJetsBdisc->push_back(dynamic_cast<const pat::Jet *>(puppiSubJets->at(ip).daughter(id))->bDiscriminator(bTagKeyString_.c_str()));
     }
  }
  //Puppi End ************
  //

  int cntJetLowPt = 0;
  for(unsigned int ij=0; ij < extJets.size(); ij++)
  {
    const pat::Jet& jet = extJets[ij];

 float trialDeepCSVb = jet.bDiscriminator((deepCSVBJetTags_+":probb").c_str());
    DeepCSVb->push_back(trialDeepCSVb);

    float trialDeepCSVc = jet.bDiscriminator((deepCSVBJetTags_+":probc").c_str());
    DeepCSVc->push_back(trialDeepCSVc);

    float trialDeepCSVl = jet.bDiscriminator((deepCSVBJetTags_+":probudsg").c_str());
    DeepCSVl->push_back(trialDeepCSVl);

    float trialDeepCSVbb = jet.bDiscriminator((deepCSVBJetTags_+":probbb").c_str());
    DeepCSVbb->push_back(trialDeepCSVbb);

    float trialDeepCSVcc = jet.bDiscriminator((deepCSVBJetTags_+":probcc").c_str());
    DeepCSVcc->push_back(trialDeepCSVcc);

    float trialDeepCSVbN = jet.bDiscriminator((deepCSVNegBJetTags_+":probb").c_str());
    DeepCSVbN->push_back(trialDeepCSVbN);

    float trialDeepCSVcN = jet.bDiscriminator((deepCSVNegBJetTags_+":probc").c_str());
    DeepCSVcN->push_back(trialDeepCSVcN);

    float trialDeepCSVlN = jet.bDiscriminator((deepCSVNegBJetTags_+":probudsg").c_str());
    DeepCSVlN->push_back(trialDeepCSVlN);

    float trialDeepCSVbbN = jet.bDiscriminator((deepCSVNegBJetTags_+":probbb").c_str());
    DeepCSVbbN->push_back(trialDeepCSVbbN);

    float trialDeepCSVccN = jet.bDiscriminator((deepCSVNegBJetTags_+":probcc").c_str());
    DeepCSVccN->push_back(trialDeepCSVccN);

    float trialDeepCSVbP = jet.bDiscriminator((deepCSVPosBJetTags_+":probb").c_str());
    DeepCSVbP->push_back(trialDeepCSVbP);

    float trialDeepCSVcP = jet.bDiscriminator((deepCSVPosBJetTags_+":probc").c_str());
DeepCSVcP->push_back(trialDeepCSVcP);

    float trialDeepCSVlP = jet.bDiscriminator((deepCSVPosBJetTags_+":probudsg").c_str());
    DeepCSVlP->push_back(trialDeepCSVlP);

    float trialDeepCSVbbP = jet.bDiscriminator((deepCSVPosBJetTags_+":probbb").c_str());
    DeepCSVbbP->push_back(trialDeepCSVbbP);

    float trialDeepCSVccP = jet.bDiscriminator((deepCSVPosBJetTags_+":probcc").c_str());
    DeepCSVccP->push_back(trialDeepCSVccP);

    float tri_CombinedIVF =jet.bDiscriminator(combinedIVFSVBJetTags_.c_str());
    CombinedIVF->push_back(tri_CombinedIVF);

    float tri_CombinedIVF_P =jet.bDiscriminator(combinedIVFSVPosBJetTags_.c_str());
    CombinedIVFP ->push_back(tri_CombinedIVF_P);

    float tri_CombinedIVF_N =jet.bDiscriminator(combinedIVFSVNegBJetTags_.c_str());
    CombinedIVFN ->push_back(tri_CombinedIVF_N);

    float Proba = jet.bDiscriminator(jetPBJetTags_.c_str());
    JetProba->push_back(Proba);

    float ProbaN = jet.bDiscriminator(jetPNegBJetTags_.c_str());
    JetProbaN->push_back(ProbaN);

    float ProbaP = jet.bDiscriminator(jetPPosBJetTags_.c_str());
    JetProbaP->push_back(ProbaP);

    float Bprob = jet.bDiscriminator(jetBPBJetTags_.c_str());
    JetBprob->push_back(Bprob);

    float BprobN = jet.bDiscriminator(jetBPNegBJetTags_.c_str());
    JetBprobN->push_back(BprobN);

    float BprobP = jet.bDiscriminator(jetBPPosBJetTags_.c_str());
    JetBprobP->push_back(BprobP);

    float Tri_CombinedSvtx = jet.bDiscriminator(combinedSVBJetTags_.c_str());
    CombinedSvtx->push_back(Tri_CombinedSvtx);

    float Tri_CombinedSvtxN = jet.bDiscriminator(combinedSVNegBJetTags_.c_str());
    CombinedSvtxN->push_back(Tri_CombinedSvtxN);

    float Tri_CombinedSvtxP = jet.bDiscriminator(combinedSVPosBJetTags_.c_str());
    CombinedSvtxP-> push_back(Tri_CombinedSvtxP);

    float tri_Svtx = jet.bDiscriminator(simpleSVHighEffBJetTags_.c_str());
    Svtx  ->push_back(tri_Svtx);
    float tri_SvtxN = jet.bDiscriminator(simpleSVNegHighEffBJetTags_.c_str());
    SvtxN ->push_back(tri_SvtxN);
    float tri_SvtxHP  = jet.bDiscriminator(simpleSVHighPurBJetTags_.c_str());
    SvtxHP->push_back(tri_SvtxHP);
    float tri_SvtxNHP = jet.bDiscriminator(simpleSVNegHighPurBJetTags_.c_str());
    SvtxNHP->push_back(tri_SvtxNHP);

    float tri_SoftM  = jet.bDiscriminator(softPFMuonBJetTags_.c_str());
    SoftM->push_back(tri_SoftM);
    float tri_SoftMN = jet.bDiscriminator(softPFMuonNegBJetTags_.c_str());
    SoftMN->push_back(tri_SoftMN);
    float tri_SoftMP = jet.bDiscriminator(softPFMuonPosBJetTags_.c_str());
    SoftMP->push_back(tri_SoftMP);

    float tri_SoftE  = jet.bDiscriminator(softPFElectronBJetTags_.c_str());
    SoftE->push_back(tri_SoftE);
    float tri_SoftEN = jet.bDiscriminator(softPFElectronNegBJetTags_.c_str());
    SoftEN->push_back(tri_SoftEN);
    float tri_SoftEP = jet.bDiscriminator(softPFElectronPosBJetTags_.c_str());
    SoftEP->push_back(tri_SoftEP);

    float tri_DoubleSV = jet.bDiscriminator(doubleSVBJetTags_.c_str());
    DoubleSV->push_back(tri_DoubleSV);

    float tri_cMVA = jet.bDiscriminator(cMVABJetTags_.c_str());
    cMVA->push_back(tri_cMVA);
    float tri_cMVAv2 = jet.bDiscriminator(cMVAv2BJetTags_.c_str());
    cMVAv2->push_back(tri_cMVAv2);
    float tri_cMVAv2Neg = jet.bDiscriminator(cMVAv2NegBJetTags_.c_str());
    cMVAv2Neg->push_back(tri_cMVAv2Neg);
    float tri_cMVAv2Pos = jet.bDiscriminator(cMVAv2PosBJetTags_.c_str());
    cMVAv2Pos->push_back(tri_cMVAv2Pos);

    float tri_CvsB = jet.bDiscriminator(CvsBCJetTags_.c_str());
    CvsB->push_back(tri_CvsB);
    float tri_CvsBNeg = jet.bDiscriminator(CvsBNegCJetTags_.c_str());
    CvsBNeg->push_back(tri_CvsBNeg);
    float tri_CvsBPos = jet.bDiscriminator(CvsBPosCJetTags_.c_str());
    CvsBPos->push_back(tri_CvsBPos);
    float tri_CvsL = jet.bDiscriminator(CvsLCJetTags_.c_str());
    CvsL->push_back(tri_CvsL);
    float tri_CvsLNeg = jet.bDiscriminator(CvsLNegCJetTags_.c_str());
    CvsLNeg->push_back(tri_CvsLNeg);
    float tri_CvsLPos = jet.bDiscriminator(CvsLPosCJetTags_.c_str());
    CvsLPos->push_back(tri_CvsLPos);

    float trialDeepFlavorb = jet.bDiscriminator((deepFlavorBJetTags_+":probb").c_str());
    DeepFlavorb->push_back(trialDeepFlavorb);

    float trialDeepFlavorbb = jet.bDiscriminator((deepFlavorBJetTags_+":probbb").c_str());
    DeepFlavorbb->push_back(trialDeepFlavorbb);

    float trialDeepFlavorlepb = jet.bDiscriminator((deepFlavorBJetTags_+":problepb").c_str());
    DeepFlavorlepb->push_back(trialDeepFlavorlepb);

    float trialDeepFlavorc = jet.bDiscriminator((deepFlavorBJetTags_+":probc").c_str());
    DeepFlavorc->push_back(trialDeepFlavorc);

    float trialDeepFlavoruds = jet.bDiscriminator((deepFlavorBJetTags_+":probuds").c_str());
    DeepFlavoruds->push_back(trialDeepFlavoruds);

    float trialDeepFlavorg = jet.bDiscriminator((deepFlavorBJetTags_+":probg").c_str());
    DeepFlavorg->push_back(trialDeepFlavorg);
    TLorentzVector perJetLVec;
    perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
    jetsLVec->push_back(perJetLVec);

    //Additional jec qualities
    std::vector<std::string> availableJECSets   = jet.availableJECSets();
    std::vector<std::string> availableJECLevels = jet.availableJECLevels(jet.currentJECSet());
    float scaleRawToFull = jet.jecFactor(jet.currentJECLevel(), "none", jet.currentJECSet())/jet.jecFactor("Uncorrected", "none", jet.currentJECSet());
    //double scaleRawToFull = jet.jecFactor(availableJECLevels.back())/jet.jecFactor("Uncorrected");
    recoJetsJecScaleRawToFull->push_back(scaleRawToFull);
    if( debug_ && ij==0 )
    {
      std::cout<<"\nAvailable JEC sets:"<<"   current : "<<jet.currentJECSet().c_str()<<std::endl;
      for(unsigned int ia=0; ia<availableJECSets.size(); ia++)
      {
         std::cout<<"ia : "<<ia<<"  --> "<<availableJECSets[ia].c_str()<<std::endl;
      }
      std::cout<<"\nAvailable JEC levels:"<<"   current : "<<jet.currentJECLevel().c_str()<<std::endl;
      for(unsigned int ia=0; ia<availableJECLevels.size(); ia++)
      {
        std::cout<<"ia : "<<ia<<"  --> "<<availableJECLevels[ia].c_str()<<std::endl;
      }
      std::cout<<"scaleRawToFull : "<<scaleRawToFull<<"  current : "<<jet.jecFactor(availableJECLevels.back())<<"  uncor : "<<jet.jecFactor("Uncorrected")<<std::endl;
    }

    //get JEC unc for this jet, using corrected pT
    jecUnc->setJetEta(jet.eta());
    jecUnc->setJetPt(jet.pt());

    float uncertainty = jecUnc->getUncertainty(true);
    //safety check if uncertainty is not available for a jet
    if( uncertainty==-999. ) uncertainty = 0;
    recoJetsJecUnc->push_back(uncertainty);

    if( perJetLVec.Pt() < jetPtCut_miniAOD_ && ij < jets->size() ) cntJetLowPt ++;

    int flavor = jet.partonFlavour();
    recoJetsFlavor->push_back(flavor);

    std::string toGetName = qgTaggerKey_+":qgLikelihood";
    if( ij >= jets->size() && qgTaggerKey_ == "QGTagger" ) toGetName = qgTaggerKey_+"Other:qgLikelihood";
    float thisqgLikelihood = jet.userFloat(toGetName.c_str());
    qgLikelihood->push_back(thisqgLikelihood);
 
    toGetName = qgTaggerKey_+":ptD";
    if( ij >= jets->size() && qgTaggerKey_ == "QGTagger" ) toGetName = qgTaggerKey_+"Other:ptD";
    float thisqgPtD = jet.userFloat(toGetName.c_str());
    qgPtD->push_back(thisqgPtD);
   
    toGetName = qgTaggerKey_+":axis2"; 
    if( ij >= jets->size() && qgTaggerKey_ == "QGTagger" ) toGetName = qgTaggerKey_+"Other:axis2";
    float thisqgAxis2 = jet.userFloat(toGetName.c_str());
    qgAxis2->push_back(thisqgAxis2);
  
    toGetName = qgTaggerKey_+":axis1"; 
    if( ij >= jets->size() && qgTaggerKey_ == "QGTagger" ) toGetName = qgTaggerKey_+"Other:axis1";
    double thisqgAxis1 = jet.userFloat(toGetName.c_str());
    qgAxis1->push_back(thisqgAxis1);

    toGetName = qgTaggerKey_+":mult"; 
    if( ij >= jets->size() && qgTaggerKey_ == "QGTagger" ) toGetName = qgTaggerKey_+"Other:mult";
    int thisqgMult = jet.userInt(toGetName.c_str());
    qgMult->push_back(thisqgMult);

    float btag = jet.bDiscriminator(bTagKeyString_.c_str());
    recoJetsBtag->push_back(btag);

    float charge = jet.jetCharge();
    recoJetsCharge->push_back(charge);

    float chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
    recoJetschargedHadronEnergyFraction->push_back( chargedHadronEnergyFraction );

    float chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
    recoJetschargedEmEnergyFraction->push_back( chargedEmEnergyFraction );

    float neutralEmEnergyFraction = jet.neutralEmEnergyFraction();
    recoJetsneutralEmEnergyFraction->push_back( neutralEmEnergyFraction );

    float muonEnergyFraction = jet.muonEnergyFraction();
    recoJetsmuonEnergyFraction->push_back( muonEnergyFraction );

    //std::cout << chargedEmEnergyFraction << std::endl;

    //const std::vector<reco::PFCandidatePtr> & constituents = jet.getPFConstituents();
    //const unsigned int numConstituents = constituents.size();
    const unsigned int numConstituents = jet.numberOfDaughters();

    for(unsigned int im=0; im < muLVec_->size(); im++)
    {
      float muEta = muLVec_->at(im).Eta(), muPhi = muLVec_->at(im).Phi();
      float mindRmuonCon = 999.;
      for (unsigned int iCon = 0; iCon < numConstituents; ++iCon)
      {
        //const reco::PFCandidatePtr& constituent = constituents[iCon];
        const reco::Candidate * constituent = jet.daughter(iCon);
        float dRmuonCon = reco::deltaR(constituent->eta(), constituent->phi(), muEta, muPhi);
        if( mindRmuonCon > dRmuonCon ){ mindRmuonCon = dRmuonCon; }
      }
      if( mindRmuonCon < deltaRcon_ ) (*muMatchedJetIdx)[im] = ij;
    }

    for(unsigned int ie=0; ie < eleLVec_->size(); ie++)
    {
      float eleEta = eleLVec_->at(ie).Eta(), elePhi = eleLVec_->at(ie).Phi();
      float mindReleCon = 999.;
      for (unsigned int iCon = 0; iCon < numConstituents; ++iCon)
      {
        //const reco::PFCandidatePtr& constituent = constituents[iCon];
        const reco::Candidate * constituent = jet.daughter(iCon);
        float dReleCon = reco::deltaR(constituent->eta(), constituent->phi(), eleEta, elePhi);
        if( mindReleCon > dReleCon ){ mindReleCon = dReleCon; }
      }
      if( mindReleCon < deltaRcon_ ) (*eleMatchedJetIdx)[ie] = ij;
    }

    for(unsigned int it=0; it < trksForIsoVetoLVec_->size(); it++)
    {
      float trkEta = trksForIsoVetoLVec_->at(it).Eta(), trkPhi = trksForIsoVetoLVec_->at(it).Phi();
      float mindRtrkCon = 999.;
      for (unsigned int iCon = 0; iCon < numConstituents; ++iCon)
      {
        //          const reco::PFCandidatePtr& constituent = constituents[iCon];
        const reco::Candidate * constituent = jet.daughter(iCon);
        float dRtrkCon = reco::deltaR(constituent->eta(), constituent->phi(), trkEta, trkPhi);
        if( mindRtrkCon > dRtrkCon ){ mindRtrkCon = dRtrkCon; }
      }
      if( mindRtrkCon < deltaRcon_ ) (*trksForIsoVetoMatchedJetIdx)[it] = ij;
    }

    for(unsigned int ist=0; ist < looseisoTrksLVec_->size(); ist++)
    {
      float isotrkEta = looseisoTrksLVec_->at(ist).Eta(), isotrkPhi = looseisoTrksLVec_->at(ist).Phi();
      float mindRisotrkCon = 999.;
      for(unsigned int iCon = 0; iCon < numConstituents; ++iCon)
      {
        //const reco::PFCandidatePtr& constituent = constituents[iCon];
        const reco::Candidate * constituent = jet.daughter(iCon);
        float dRisotrkCon = reco::deltaR(constituent->eta(), constituent->phi(), isotrkEta, isotrkPhi);
        if( mindRisotrkCon > dRisotrkCon ){ mindRisotrkCon = dRisotrkCon; }
      }
      if( mindRisotrkCon < deltaRcon_ ) (*looseisoTrksMatchedJetIdx)[ist] = ij;
    }
  }

  if( cntJetLowPt ) std::cout<<"WARNING ... NON ZERO ("<<cntJetLowPt<<") number of jets with pt < "<<jetPtCut_miniAOD_<<std::endl;

  std::unique_ptr<int> nJets (new int);

  *nJets = jetsLVec->size();

  // store in the event
  // iEvent.put(prod);
  iEvent.put(std::move(jetsLVec), "jetsLVec");
  iEvent.put(std::move(recoJetsFlavor), "recoJetsFlavor");
  iEvent.put(std::move(recoJetsBtag), "recoJetsBtag");
  iEvent.put(std::move(recoJetsCharge), "recoJetsCharge");
  iEvent.put(std::move(recoJetsJecUnc), "recoJetsJecUnc");
  iEvent.put(std::move(recoJetsJecScaleRawToFull), "recoJetsJecScaleRawToFull");
  iEvent.put(std::move(nJets), "nJets");

  iEvent.put(std::move(qgLikelihood), "qgLikelihood");
  iEvent.put(std::move(qgPtD), "qgPtD");
  iEvent.put(std::move(qgAxis2), "qgAxis2");
  iEvent.put(std::move(qgAxis1), "qgAxis1");
  iEvent.put(std::move(qgMult), "qgMult");

  iEvent.put(std::move(recoJetschargedHadronEnergyFraction), "recoJetschargedHadronEnergyFraction");
  iEvent.put(std::move(recoJetschargedEmEnergyFraction), "recoJetschargedEmEnergyFraction");
  iEvent.put(std::move(recoJetsneutralEmEnergyFraction), "recoJetsneutralEmEnergyFraction");

  iEvent.put(std::move(recoJetsmuonEnergyFraction), "recoJetsmuonEnergyFraction");

  iEvent.put(std::move(muMatchedJetIdx), "muMatchedJetIdx");
  iEvent.put(std::move(eleMatchedJetIdx), "eleMatchedJetIdx");

  iEvent.put(std::move(trksForIsoVetoMatchedJetIdx), "trksForIsoVetoMatchedJetIdx");
  iEvent.put(std::move(looseisoTrksMatchedJetIdx), "looseisoTrksMatchedJetIdx");
  //PUPPI
  iEvent.put(std::move(puppiJetsLVec), "puppiJetsLVec");
  iEvent.put(std::move(puppiSubJetsLVec), "puppiSubJetsLVec");
  iEvent.put(std::move(puppiSubJetsBdisc), "puppiSubJetsBdisc");
  iEvent.put(std::move(puppisoftDropMass),"puppisoftDropMass");
  iEvent.put(std::move(puppitau1),"puppitau1");
  iEvent.put(std::move(puppitau2),"puppitau2");
  iEvent.put(std::move(puppitau3),"puppitau3");

  iEvent.put(std::move(DeepCSVb),"DeepCSVb");
  iEvent.put(std::move(DeepCSVc),"DeepCSVc");
  iEvent.put(std::move(DeepCSVl),"DeepCSVl");
  iEvent.put(std::move(DeepCSVbb),"DeepCSVbb");
  iEvent.put(std::move(DeepCSVcc),"DeepCSVcc");

  iEvent.put(std::move(DeepCSVbN),"DeepCSVbN");
  iEvent.put(std::move(DeepCSVcN),"DeepCSVcN");
  iEvent.put(std::move(DeepCSVlN),"DeepCSVlN");
  iEvent.put(std::move(DeepCSVbbN),"DeepCSVbbN");
  iEvent.put(std::move(DeepCSVccN),"DeepCSVccN");

  iEvent.put(std::move(DeepCSVbP),"DeepCSVbP");
  iEvent.put(std::move(DeepCSVcP),"DeepCSVcP");
  iEvent.put(std::move(DeepCSVlP),"DeepCSVlP");
  iEvent.put(std::move(DeepCSVbbP),"DeepCSVbbP");
  iEvent.put(std::move(DeepCSVccP),"DeepCSVccP");

  iEvent.put(std::move(DeepFlavorb), "DeepFlavorb");
  iEvent.put(std::move(DeepFlavorbb), "DeepFlavorbb");
  iEvent.put(std::move(DeepFlavorlepb), "DeepFlavorlepb");
  iEvent.put(std::move(DeepFlavorc), "DeepFlavorc");
  iEvent.put(std::move(DeepFlavoruds), "DeepFlavoruds");
  iEvent.put(std::move(DeepFlavorg), "DeepFlavorg");

  iEvent.put(std::move(CombinedSvtx),"CombinedSvtx");
  iEvent.put(std::move(CombinedSvtxN),"CombinedSvtxN");
  iEvent.put(std::move(CombinedSvtxP),"CombinedSvtxP");

  iEvent.put(std::move(CombinedIVF),"CombinedIVF");
  iEvent.put(std::move(CombinedIVFP),"CombinedIVFP");
  iEvent.put(std::move(CombinedIVFN),"CombinedIVFN");

  iEvent.put(std::move(Svtx),"Svtx");
  iEvent.put(std::move(SvtxN),"SvtxN");
  iEvent.put(std::move(SvtxHP),"SvtxHP");
  iEvent.put(std::move(SvtxNHP),"SvtxNHP");

  iEvent.put(std::move(SoftM),"SoftM");
  iEvent.put(std::move(SoftMN),"SoftMN");
  iEvent.put(std::move(SoftMP),"SoftMP");

  iEvent.put(std::move(SoftE),"SoftE");
  iEvent.put(std::move(SoftEN),"SoftEN");
  iEvent.put(std::move(SoftEP),"SoftEP");

  iEvent.put(std::move(DoubleSV),"DoubleSV");

  iEvent.put(std::move(cMVA),"cMVA");
  iEvent.put(std::move(cMVAv2),"cMVAv2");
  iEvent.put(std::move(cMVAv2Neg),"cMVAv2Neg");
  iEvent.put(std::move(cMVAv2Pos),"cMVAv2Pos");

  iEvent.put(std::move(CvsB),"CvsB");
  iEvent.put(std::move(CvsBNeg),"CvsBNeg");
  iEvent.put(std::move(CvsBPos),"CvsBPos");
  iEvent.put(std::move(CvsL),"CvsL");
  iEvent.put(std::move(CvsLNeg),"CvsLNeg");
  iEvent.put(std::move(CvsLPos),"CvsLPos");

  return true;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(prodJets);
