
#include <memory>
#include <algorithm>
#include <vector>
#include <map>

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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

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

  edm::InputTag jetSrc_;
  std::string bTagKeyString_;
  bool debug_;

  double jetPtCut_miniAOD_, genMatch_dR_;

  edm::EDGetTokenT<std::vector<pat::Jet> >JetTok_;

  edm::EDGetTokenT<std::vector<TLorentzVector> >GenDecayLVec_Tok_;
  edm::EDGetTokenT<std::vector<int> >GenDecayMomRefVec_Tok_;

  edm::EDGetTokenT<std::vector<TLorentzVector> >EleLVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >MuLVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >TrksForIsoVetolVec_Tok_;
  edm::EDGetTokenT<std::vector<TLorentzVector> >LooseIsoTrksVec_Tok_;

  edm::EDGetTokenT<std::vector<pat::Jet>> PuppiJetsSrc_Tok_;
  edm::EDGetTokenT<std::vector<pat::Jet>> PuppiSubJetsSrc_Tok_;

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

  std::string jetBPBJetTags_;

  std::string deepCSVBJetTags_;

  std::string deepFlavorBJetTags_;

  std::string   CvsBCJetTags_;
  std::string   CvsLCJetTags_;
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
                     if( (p->dxy()*p->dxy()) / (p->dxyError()*p->dxyError()) < 25. ) ++totalMult_; // this cut only      applies to multiplicity
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
  jetSrc_      = iConfig.getParameter<edm::InputTag>("jetSrc");
  bTagKeyString_ = iConfig.getParameter<std::string>("bTagKeyString");

  jetPBJetTags_        = iConfig.getParameter<std::string>("jetPBJetTags");

  jetBPBJetTags_        = iConfig.getParameter<std::string>("jetBPBJetTags");

  deepCSVBJetTags_    = iConfig.getParameter<std::string>("deepCSVBJetTags");

  deepFlavorBJetTags_    = iConfig.getParameter<std::string>("deepFlavorBJetTags");

  CvsBCJetTags_             = iConfig.getParameter<std::string>("CvsBCJetTags");
  CvsLCJetTags_             = iConfig.getParameter<std::string>("CvsLCJetTags");

  debug_       = iConfig.getParameter<bool>("debug");

  jetPtCut_miniAOD_ = iConfig.getUntrackedParameter<double>("jetPtCut_miniAOD", 10);
  genMatch_dR_ = iConfig.getUntrackedParameter<double>("genMatch_dR", 1.0);

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
  EleLVec_Tok_=consumes<std::vector<TLorentzVector> >(eleLVec_Src_);
  MuLVec_Tok_=consumes<std::vector<TLorentzVector> >(muLVec_Src_);
  TrksForIsoVetolVec_Tok_=consumes<std::vector<TLorentzVector> >(trksForIsoVetoLVec_Src_);
  LooseIsoTrksVec_Tok_=consumes<std::vector<TLorentzVector> >(looseisoTrksLVec_Src_);
  PuppiJetsSrc_Tok_ = consumes<std::vector<pat::Jet>>(puppiJetsSrc_);
  PuppiSubJetsSrc_Tok_ = consumes<std::vector<pat::Jet>>(puppiSubJetsSrc_);

  //produces<std::vector<pat::Jet> >("");
  produces<std::vector<TLorentzVector> >("jetsLVec");
  produces<std::vector<int> >("recoJetsFlavor");
  produces<std::vector<float> >("recoJetsCSVv2");
  produces<std::vector<float> >("recoJetsCharge");
  produces<std::vector<float> >("recoJetsJecUnc");
  produces<std::vector<float> >("recoJetsJecScaleRawToFull");
  produces<int>("nJets");
  produces<std::vector<float> >("qgLikelihood");
  produces<std::vector<float> >("qgPtD");
  produces<std::vector<float> >("qgAxis2");
  produces<std::vector<float> >("qgAxis1");
  produces<std::vector<float> >("qgMult");

  //produce variables needed for Lost Lepton study, added by hua.wei@cern.ch
  produces<std::vector<float> >("recoJetschargedHadronEnergyFraction");
  produces<std::vector<float> >("recoJetschargedEmEnergyFraction");
  produces<std::vector<float> >("recoJetsneutralEmEnergyFraction");
  produces<std::vector<float> >("recoJetsHFHadronEnergyFraction");
  produces<std::vector<float> >("recoJetsmuonEnergyFraction");
  produces<std::vector<float> >("recoJetsneutralEnergyFraction");
  produces<std::vector<float> >("recoJetsHFEMEnergyFraction");
  
  produces<std::vector<float> >("PhotonEnergyFraction");
  produces<std::vector<float> >("ElectronEnergyFraction");

  produces<std::vector<int> >("muMatchedJetIdx");
  produces<std::vector<int> >("eleMatchedJetIdx");

  produces<std::vector<float> >("ChargedHadronMultiplicity");
  produces<std::vector<float> >("NeutralHadronMultiplicity");
  produces<std::vector<float> >("PhotonMultiplicity");
  produces<std::vector<float> >("ElectronMultiplicity");
  produces<std::vector<float> >("MuonMultiplicity");

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
  produces<std::vector<float> >("JetBprob");

  produces<std::vector<float> >("CombinedSvtx");

  produces<std::vector<float> >("DeepCSVb");
  produces<std::vector<float> >("DeepCSVc");
  produces<std::vector<float> >("DeepCSVl");
  produces<std::vector<float> >("DeepCSVbb");
  produces<std::vector<float> >("DeepCSVcc");

  produces<std::vector<float> >("DeepFlavorb");
  produces<std::vector<float> >("DeepFlavorbb");
  produces<std::vector<float> >("DeepFlavorlepb");
  produces<std::vector<float> >("DeepFlavorc");
  produces<std::vector<float> >("DeepFlavoruds");
  produces<std::vector<float> >("DeepFlavorg");

  produces<std::vector<float> >("CversusB");
  produces<std::vector<float> >("CversusL");

}


prodJets::~prodJets() 
{
}


bool prodJets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  edm::Handle<std::vector<pat::Jet> > jets; 
  iEvent.getByToken(JetTok_, jets);

  //get the JEC uncertainties
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jetType_, JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  std::unique_ptr<JetCorrectionUncertainty> jecUnc( new JetCorrectionUncertainty(JetCorPar) );

  iEvent.getByToken(EleLVec_Tok_, eleLVec_);
  iEvent.getByToken(MuLVec_Tok_, muLVec_);

  iEvent.getByToken(TrksForIsoVetolVec_Tok_, trksForIsoVetoLVec_);
  iEvent.getByToken(LooseIsoTrksVec_Tok_,looseisoTrksLVec_);

  std::vector<pat::Jet> extJets = (*jets);
  //PUPPI
  iEvent.getByToken(PuppiJetsSrc_Tok_, puppiJets);
  iEvent.getByToken(PuppiSubJetsSrc_Tok_, puppiSubJets);

  //check which ones to keep
  std::unique_ptr<std::vector<TLorentzVector> > jetsLVec(new std::vector<TLorentzVector>());
  std::unique_ptr<std::vector<int> > recoJetsFlavor(new std::vector<int>());
  std::unique_ptr<std::vector<float> > recoJetsCSVv2(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsCharge(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsJecUnc(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsJecScaleRawToFull(new std::vector<float>());

  std::unique_ptr<std::vector<float> > JetProba(new std::vector<float>());
  std::unique_ptr<std::vector<float> > JetBprob(new std::vector<float>());

  std::unique_ptr<std::vector<float> > DeepCSVb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVc(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVl(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVbb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepCSVcc(new std::vector<float>());

  std::unique_ptr<std::vector<float> > DeepFlavorb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorbb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorlepb(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorc(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavoruds(new std::vector<float>());
  std::unique_ptr<std::vector<float> > DeepFlavorg(new std::vector<float>());

  std::unique_ptr<std::vector<float> >CversusB(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CversusL(new std::vector<float>());

  std::unique_ptr<std::vector<float> > qgLikelihood(new std::vector<float>());
  std::unique_ptr<std::vector<float> > qgPtD(new std::vector<float>());
  std::unique_ptr<std::vector<float> > qgAxis2(new std::vector<float>());
  std::unique_ptr<std::vector<float> > qgAxis1(new std::vector<float>());
  std::unique_ptr<std::vector<float> > qgMult(new std::vector<float>());

  std::unique_ptr<std::vector<float> > recoJetschargedHadronEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetschargedEmEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsneutralEmEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsmuonEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsneutralEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsHFEMEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > recoJetsHFHadronEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > PhotonEnergyFraction(new std::vector<float>());
  std::unique_ptr<std::vector<float> > ElectronEnergyFraction(new std::vector<float>());

  std::unique_ptr<std::vector<float> > ChargedHadronMultiplicity(new std::vector<float>());
  std::unique_ptr<std::vector<float> > NeutralHadronMultiplicity(new std::vector<float>());
  std::unique_ptr<std::vector<float> > PhotonMultiplicity(new std::vector<float>());
  std::unique_ptr<std::vector<float> > ElectronMultiplicity(new std::vector<float>());
  std::unique_ptr<std::vector<float> > MuonMultiplicity(new std::vector<float>());

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
    //std::cout<<trialDeepCSVb<<std::endl;
    DeepCSVb->push_back(trialDeepCSVb);

    float trialDeepCSVc = jet.bDiscriminator((deepCSVBJetTags_+":probc").c_str());
    DeepCSVc->push_back(trialDeepCSVc);

    float trialDeepCSVl = jet.bDiscriminator((deepCSVBJetTags_+":probudsg").c_str());
    DeepCSVl->push_back(trialDeepCSVl);

    float trialDeepCSVbb = jet.bDiscriminator((deepCSVBJetTags_+":probbb").c_str());
    DeepCSVbb->push_back(trialDeepCSVbb);

    float trialDeepCSVcc = jet.bDiscriminator((deepCSVBJetTags_+":probcc").c_str());
    DeepCSVcc->push_back(trialDeepCSVcc);

    float Proba = jet.bDiscriminator(jetPBJetTags_.c_str());
    JetProba->push_back(Proba);

    float Bprob = jet.bDiscriminator(jetBPBJetTags_.c_str());
    JetBprob->push_back(Bprob);

    float tri_CvsB = jet.bDiscriminator(CvsBCJetTags_.c_str());
    CversusB->push_back(tri_CvsB);
    float tri_CvsL = jet.bDiscriminator(CvsLCJetTags_.c_str());
    CversusL->push_back(tri_CvsL);

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
    float scaleRawToFull = jet.jecFactor(jet.currentJECLevel(), "none", jet.currentJECSet())/jet.jecFactor("Uncorrected", "none", jet.currentJECSet());
    //double scaleRawToFull = jet.jecFactor(availableJECLevels.back())/jet.jecFactor("Uncorrected");
    recoJetsJecScaleRawToFull->push_back(scaleRawToFull);
    if( debug_ && ij==0 )
    {
      std::vector<std::string> availableJECSets   = jet.availableJECSets();
      std::vector<std::string> availableJECLevels = jet.availableJECLevels(jet.currentJECSet());
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
    float thisqgAxis1 = jet.userFloat(toGetName.c_str());
    //std::cout<<thisqgAxis1<<std::endl;
    qgAxis1->push_back(thisqgAxis1);
    //std::cout<<qgAxis1->at(0)<<std::endl;

    toGetName = qgTaggerKey_+":mult"; 
    if( ij >= jets->size() && qgTaggerKey_ == "QGTagger" ) toGetName = qgTaggerKey_+"Other:mult";
    int thisqgMult = jet.userInt(toGetName.c_str());
    qgMult->push_back(thisqgMult);

    float btag = jet.bDiscriminator(bTagKeyString_.c_str());
    recoJetsCSVv2->push_back(btag);

    float charge = jet.jetCharge();
    recoJetsCharge->push_back(charge);

    float chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
    recoJetschargedHadronEnergyFraction->push_back( chargedHadronEnergyFraction );

    double neutralHadronEnergyFraction = jet.neutralHadronEnergyFraction();
    recoJetsneutralEnergyFraction->push_back( neutralHadronEnergyFraction );

    float photonEnergyFraction = jet.photonEnergyFraction();
    PhotonEnergyFraction->push_back( photonEnergyFraction );
    
    float electronEnergyFraction = jet.electronEnergyFraction();
    ElectronEnergyFraction->push_back( electronEnergyFraction );
    
    recoJetsHFHadronEnergyFraction->push_back(jet.HFHadronEnergyFraction());
    recoJetsHFEMEnergyFraction->push_back(jet.HFEMEnergyFraction());
    
    float chargedHadronMultiplicity = jet.chargedHadronMultiplicity();
    ChargedHadronMultiplicity->push_back( chargedHadronMultiplicity );
    
    float neutralHadronMultiplicity = jet.neutralHadronMultiplicity();
    NeutralHadronMultiplicity->push_back( neutralHadronMultiplicity );
    
    float photonMultiplicity = jet.photonMultiplicity();
    PhotonMultiplicity->push_back( photonMultiplicity );
    
    float electronMultiplicity = jet.electronMultiplicity();
    ElectronMultiplicity->push_back( electronMultiplicity );

    double muonMultiplicity1 = jet.muonMultiplicity();
    MuonMultiplicity->push_back( muonMultiplicity1 );

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
  iEvent.put(std::move(recoJetsCSVv2), "recoJetsCSVv2");
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
  
  iEvent.put(std::move(ChargedHadronMultiplicity), "ChargedHadronMultiplicity");
  iEvent.put(std::move(NeutralHadronMultiplicity), "NeutralHadronMultiplicity");
  iEvent.put(std::move(PhotonMultiplicity), "PhotonMultiplicity");
  iEvent.put(std::move(ElectronMultiplicity), "ElectronMultiplicity");
  iEvent.put(std::move(MuonMultiplicity), "MuonMultiplicity");
  
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

  iEvent.put(std::move(recoJetsneutralEnergyFraction), "recoJetsneutralEnergyFraction");
  iEvent.put(std::move(recoJetsHFHadronEnergyFraction), "recoJetsHFHadronEnergyFraction");
  iEvent.put(std::move(recoJetsHFEMEnergyFraction), "recoJetsHFEMEnergyFraction");

  iEvent.put(std::move(PhotonEnergyFraction), "PhotonEnergyFraction");
  iEvent.put(std::move(ElectronEnergyFraction), "ElectronEnergyFraction");
  
  iEvent.put(std::move(DeepFlavorb), "DeepFlavorb");
  iEvent.put(std::move(DeepFlavorbb), "DeepFlavorbb");
  iEvent.put(std::move(DeepFlavorlepb), "DeepFlavorlepb");
  iEvent.put(std::move(DeepFlavorc), "DeepFlavorc");
  iEvent.put(std::move(DeepFlavoruds), "DeepFlavoruds");
  iEvent.put(std::move(DeepFlavorg), "DeepFlavorg");

  iEvent.put(std::move(CversusB),"CversusB");
  iEvent.put(std::move(CversusL),"CversusL");

  return true;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(prodJets);
