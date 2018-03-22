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

//#include "NNKit/FatJetNN/interface/FatJetNN.h"
//
class prodBTag : public edm::EDFilter 
{
public:

    explicit prodBTag(const edm::ParameterSet & iConfig);
    ~prodBTag();


private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);
  edm::EDGetTokenT<std::vector<pat::Jet> >JetTok_;
 
  bool isData_;
  
  edm::InputTag jetSrc_, jetOtherSrc_;
  edm::Handle<std::vector<pat::Jet> > jets, otherjets; 

  //std::string bTagKeyString_;

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

prodBTag::prodBTag(const edm::ParameterSet & iConfig) 
 {
  isData_ = true;

  //bTagKeyString_ = iConfig.getParameter<std::string>("bTagKeyString");
  jetSrc_      = iConfig.getParameter<edm::InputTag>("jetSrc");

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

  JetTok_ = consumes<std::vector<pat::Jet> >(jetSrc_);

produces<std::vector<double> >("recoJetsBtag");
produces<std::vector<double> >("JetProba");
  produces<std::vector<double> >("JetProbaN");
  produces<std::vector<double> >("JetProbaP");
  produces<std::vector<double> >("JetBprob");
  produces<std::vector<double> >("JetBprobP");
  produces<std::vector<double> >("JetBprobN");

  produces<std::vector<double> >("CombinedSvtx");
  produces<std::vector<double> >("CombinedSvtxN");
  produces<std::vector<double> >("CombinedSvtxP");

produces<std::vector<double> >("DeepCSVb");
  produces<std::vector<double> >("DeepCSVc");
  produces<std::vector<double> >("DeepCSVl");
  produces<std::vector<double> >("DeepCSVbb");
  produces<std::vector<double> >("DeepCSVcc");

  produces<std::vector<double> >("DeepCSVbN");
  produces<std::vector<double> >("DeepCSVcN");
  produces<std::vector<double> >("DeepCSVlN");
  produces<std::vector<double> >("DeepCSVbbN");
  produces<std::vector<double> >("DeepCSVccN");

  produces<std::vector<double> >("DeepCSVbP");
  produces<std::vector<double> >("DeepCSVcP");
  produces<std::vector<double> >("DeepCSVlP");
  produces<std::vector<double> >("DeepCSVbbP");
  produces<std::vector<double> >("DeepCSVccP");

  produces<std::vector<double> >("DeepFlavorb");
  produces<std::vector<double> >("DeepFlavorbb");
  produces<std::vector<double> >("DeepFlavorlepb");
  produces<std::vector<double> >("DeepFlavorc");
  produces<std::vector<double> >("DeepFlavoruds");
  produces<std::vector<double> >("DeepFlavorg");

  produces<std::vector<float> >("CombinedIVF");
  produces<std::vector<float> >("CombinedIVFP");
  produces<std::vector<float> >("CombinedIVFN");

  produces<std::vector<double> >("Svtx");
  produces<std::vector<double> >("SvtxN");
  produces<std::vector<double> >("SvtxHP");
  produces<std::vector<double> >("SvtxNHP");

  produces<std::vector<double> >("SoftM");
  produces<std::vector<double> >("SoftMN");
  produces<std::vector<double> >("SoftMP");

  produces<std::vector<double> >("SoftE");
  produces<std::vector<double> >("SoftEN");
  produces<std::vector<double> >("SoftEP");

  produces<std::vector<double> >("DoubleSV");

  produces<std::vector<double> >("cMVA");
  produces<std::vector<double> >("cMVAv2");
  produces<std::vector<double> >("cMVAv2Neg");
  produces<std::vector<double> >("cMVAv2Pos");

  produces<std::vector<double> >("CvsB");
  produces<std::vector<double> >("CvsBNeg");
  produces<std::vector<double> >("CvsBPos");
  produces<std::vector<double> >("CvsL");
  produces<std::vector<double> >("CvsLNeg");
  produces<std::vector<double> >("CvsLPos");

}


prodBTag::~prodBTag() 
{
}

bool prodBTag::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
 if( !iEvent.isRealData() ) isData_ = false;

 iEvent.getByToken(JetTok_, jets);
std::vector<pat::Jet> extJets = (*jets);

  std::auto_ptr<std::vector<double> > recoJetsBtag(new std::vector<double>());

  std::unique_ptr<std::vector<double> > JetProba(new std::vector<double>());
  std::unique_ptr<std::vector<double> > JetProbaN(new std::vector<double>());
  std::unique_ptr<std::vector<double> > JetProbaP(new std::vector<double>());
  std::unique_ptr<std::vector<double> > JetBprob(new std::vector<double>());
  std::unique_ptr<std::vector<double> > JetBprobP(new std::vector<double>());
  std::unique_ptr<std::vector<double> > JetBprobN(new std::vector<double>());

  std::unique_ptr<std::vector<double> >CombinedSvtx(new std::vector<double>());
  std::unique_ptr<std::vector<double> >CombinedSvtxN(new std::vector<double>());
  std::unique_ptr<std::vector<double> >CombinedSvtxP(new std::vector<double>());

  std::unique_ptr<std::vector<double> > DeepCSVb(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVc(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVl(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVbb(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVcc(new std::vector<double>());

  std::unique_ptr<std::vector<double> > DeepCSVbN(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVcN(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVlN(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVbbN(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVccN(new std::vector<double>());

  std::unique_ptr<std::vector<double> > DeepCSVbP(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVcP(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVlP(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVbbP(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepCSVccP(new std::vector<double>());

  std::unique_ptr<std::vector<double> > DeepFlavorb(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepFlavorbb(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepFlavorlepb(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepFlavorc(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepFlavoruds(new std::vector<double>());
  std::unique_ptr<std::vector<double> > DeepFlavorg(new std::vector<double>());

  std::unique_ptr<std::vector<float> >CombinedIVF(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CombinedIVFN(new std::vector<float>());
  std::unique_ptr<std::vector<float> >CombinedIVFP(new std::vector<float>());

  std::unique_ptr<std::vector<double> >Svtx(new std::vector<double>());
  std::unique_ptr<std::vector<double> >SvtxN(new std::vector<double>());
  std::unique_ptr<std::vector<double> >SvtxHP(new std::vector<double>());
  std::unique_ptr<std::vector<double> >SvtxNHP(new std::vector<double>());

  std::unique_ptr<std::vector<double> >SoftM(new std::vector<double>());
  std::unique_ptr<std::vector<double> >SoftMN(new std::vector<double>());
  std::unique_ptr<std::vector<double> >SoftMP(new std::vector<double>());

  std::unique_ptr<std::vector<double> >SoftE(new std::vector<double>());
  std::unique_ptr<std::vector<double> >SoftEN(new std::vector<double>());
  std::unique_ptr<std::vector<double> >SoftEP(new std::vector<double>());

  std::unique_ptr<std::vector<double> >DoubleSV(new std::vector<double>());

  std::unique_ptr<std::vector<double> >cMVA(new std::vector<double>());
  std::unique_ptr<std::vector<double> >cMVAv2(new std::vector<double>());
  std::unique_ptr<std::vector<double> >cMVAv2Neg(new std::vector<double>());
  std::unique_ptr<std::vector<double> >cMVAv2Pos(new std::vector<double>());

  std::unique_ptr<std::vector<double> >CvsB(new std::vector<double>());
  std::unique_ptr<std::vector<double> >CvsBNeg(new std::vector<double>());
  std::unique_ptr<std::vector<double> >CvsBPos(new std::vector<double>());
  std::unique_ptr<std::vector<double> >CvsL(new std::vector<double>());
  std::unique_ptr<std::vector<double> >CvsLNeg(new std::vector<double>());
  std::unique_ptr<std::vector<double> >CvsLPos(new std::vector<double>());

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

    double trialDeepCSVccN = jet.bDiscriminator((deepCSVNegBJetTags_+":probcc").c_str());
    DeepCSVccN->push_back(trialDeepCSVccN);

    double trialDeepCSVbP = jet.bDiscriminator((deepCSVPosBJetTags_+":probb").c_str());
    DeepCSVbP->push_back(trialDeepCSVbP);

    double trialDeepCSVcP = jet.bDiscriminator((deepCSVPosBJetTags_+":probc").c_str());
    DeepCSVcP->push_back(trialDeepCSVcP);

    double trialDeepCSVlP = jet.bDiscriminator((deepCSVPosBJetTags_+":probudsg").c_str());
    DeepCSVlP->push_back(trialDeepCSVlP);

    double trialDeepCSVbbP = jet.bDiscriminator((deepCSVPosBJetTags_+":probbb").c_str());
    DeepCSVbbP->push_back(trialDeepCSVbbP);

    double trialDeepCSVccP = jet.bDiscriminator((deepCSVPosBJetTags_+":probcc").c_str());
    DeepCSVccP->push_back(trialDeepCSVccP);

    float tri_CombinedIVF =jet.bDiscriminator(combinedIVFSVBJetTags_.c_str());
    CombinedIVF->push_back(tri_CombinedIVF);

    float tri_CombinedIVF_P =jet.bDiscriminator(combinedIVFSVPosBJetTags_.c_str());
    CombinedIVFP ->push_back(tri_CombinedIVF_P);

    float tri_CombinedIVF_N =jet.bDiscriminator(combinedIVFSVNegBJetTags_.c_str());
    CombinedIVFN ->push_back(tri_CombinedIVF_N);

    //double btag = jet.bDiscriminator(bTagKeyString_.c_str());
    //recoJetsBtag->push_back(btag);

    double Proba = jet.bDiscriminator(jetPBJetTags_.c_str());
    JetProba->push_back(Proba);

    double ProbaN = jet.bDiscriminator(jetPNegBJetTags_.c_str());
    JetProbaN->push_back(ProbaN);

    double ProbaP = jet.bDiscriminator(jetPPosBJetTags_.c_str());
    JetProbaP->push_back(ProbaP);

    double Bprob = jet.bDiscriminator(jetBPBJetTags_.c_str());
    JetBprob->push_back(Bprob);

    double BprobN = jet.bDiscriminator(jetBPNegBJetTags_.c_str());
    JetBprobN->push_back(BprobN);

    double BprobP = jet.bDiscriminator(jetBPPosBJetTags_.c_str());
    JetBprobP->push_back(BprobP);

    double Tri_CombinedSvtx = jet.bDiscriminator(combinedSVBJetTags_.c_str());
    CombinedSvtx->push_back(Tri_CombinedSvtx);

    double Tri_CombinedSvtxN = jet.bDiscriminator(combinedSVNegBJetTags_.c_str());
    CombinedSvtxN->push_back(Tri_CombinedSvtxN);
    
    double Tri_CombinedSvtxP = jet.bDiscriminator(combinedSVPosBJetTags_.c_str());
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
}

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

DEFINE_FWK_MODULE(prodBTag);
