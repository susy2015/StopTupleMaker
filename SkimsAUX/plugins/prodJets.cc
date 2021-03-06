#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "TLorentzVector.h"

#include "NNKit/FatJetNN/interface/FatJetNN.h"
#include "NNKit/FatJetNN/interface/FatJetNNHelper.h"

class prodJets : public edm::EDFilter 
{
public:

    explicit prodJets(const edm::ParameterSet & iConfig);
    ~prodJets();

private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

    void produceAK4JetVariables(edm::Event & iEvent, const edm::EventSetup & iSetup);
    void produceAK8JetVariables(edm::Event & iEvent, const edm::EventSetup & iSetup);

    void computeJetShape(const reco::Jet * jet, bool isReco, double& totalMult_, double& ptD_, double& axis1_, double& axis2_);

    template<typename T>
    void drMatchObject(const T& objLVec, const pat::Jet& jet, unsigned int ij, std::unique_ptr<std::vector<int> >& matchVec)
    {
        const unsigned int numConstituents = jet.numberOfDaughters();
        for(unsigned int im=0; im < objLVec->size(); im++)
        {
            float muEta = (*objLVec)[im].Eta(), muPhi = (*objLVec)[im].Phi();
            float mindRobjCon = 999.;
            for (unsigned int iCon = 0; iCon < numConstituents; ++iCon)
            {
                const reco::Candidate * constituent = jet.daughter(iCon);
                float dRobjCon = reco::deltaR(constituent->eta(), constituent->phi(), muEta, muPhi);
                if( mindRobjCon > dRobjCon ){ mindRobjCon = dRobjCon; }
            }
            if( mindRobjCon < deltaRcon_ ) (*matchVec)[im] = ij;
        }
    }

    //Member variables  
    bool debug_;

    double genMatch_dR_;

    double deltaRcon_;

    std::string qgTaggerKey_;

    std::string bTagKeyString_;
    std::string deepCSVBJetTags_;
    std::string deepFlavorBJetTags_;

    std::string CvsBCJetTags_;
    std::string CvsLCJetTags_;

    std::string jetType_;

    std::string NjettinessAK8Puppi_label_;
    std::string ak8PFJetsPuppi_label_;

    std::string jetPBJetTags_;
    std::string jetBPBJetTags_;

    edm::EDGetTokenT<std::vector<pat::Jet> > JetTok_;

    edm::EDGetTokenT<std::vector<TLorentzVector> > EleLVec_Tok_;
    edm::EDGetTokenT<std::vector<TLorentzVector> > MuLVec_Tok_;
    edm::EDGetTokenT<std::vector<TLorentzVector> > TrksForIsoVetolVec_Tok_;
    edm::EDGetTokenT<std::vector<TLorentzVector> > LooseIsoTrksVec_Tok_;

    edm::EDGetTokenT<std::vector<pat::Jet> > PuppiJetsSrc_Tok_;

    float AK4minPTcut_;
    float AK8minPTcut_;

    std::unique_ptr<deepntuples::FatJetNN> fatjetNN_;
};

void prodJets::computeJetShape(const reco::Jet * jet, bool isReco, double& totalMult, double& ptD, double& axis1, double& axis2)
{
    totalMult = 0;
    ptD       = 0;
    axis1     = 0;
    axis2     = 0;

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
                    if( (p->dxy()*p->dxy()) / (p->dxyError()*p->dxyError()) < 25. ) ++totalMult; // this cut only      applies to multiplicity
                } else ++totalMult;
            } else ++totalMult;
        } else { // neutral particles
            if(part->pt() < 1.0) continue;
            ++totalMult;
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
        ptD = sqrt(sum_weight)/sum_pt;
        ave_dEta  = sum_dEta  / sum_weight;
        ave_dPhi  = sum_dPhi  / sum_weight;
        ave_dEta2 = sum_dEta2 / sum_weight;
        ave_dPhi2 = sum_dPhi2 / sum_weight;
        a = ave_dEta2 - ave_dEta*ave_dEta;
        b = ave_dPhi2 - ave_dPhi*ave_dPhi;
        c = -(sum_dEta_dPhi/sum_weight - ave_dEta*ave_dPhi);
    } else ptD = 0;

    float delta = sqrt(fabs( (a-b)*(a-b) + 4*c*c ));
    if(a+b-delta > 0) axis2 = sqrt(0.5*(a+b-delta));
    else              axis2 = 0.0;
    if(a+b+delta > 0) axis1 = sqrt(0.5*(a+b+delta));
    else              axis1 = 0.0;
}

prodJets::prodJets(const edm::ParameterSet & iConfig) 
{
    //Global parameters 
    debug_                    = iConfig.getParameter<bool>("debug");

    genMatch_dR_              = iConfig.getUntrackedParameter<double>("genMatch_dR", 1.0);

    deltaRcon_                = iConfig.getUntrackedParameter<double>("deltaRcon", 0.01);

    jetType_                  = iConfig.getParameter<std::string>("jetType");

    qgTaggerKey_              = iConfig.getParameter<std::string>("qgTaggerKey");
    
    bTagKeyString_            = iConfig.getParameter<std::string>("bTagKeyString");
    deepCSVBJetTags_          = iConfig.getParameter<std::string>("deepCSVBJetTags");
    deepFlavorBJetTags_       = iConfig.getParameter<std::string>("deepFlavorBJetTags");

    CvsBCJetTags_             = iConfig.getParameter<std::string>("CvsBCJetTags");
    CvsLCJetTags_             = iConfig.getParameter<std::string>("CvsLCJetTags");

    jetPBJetTags_             = iConfig.getParameter<std::string>("jetPBJetTags");
    jetBPBJetTags_            = iConfig.getParameter<std::string>("jetBPBJetTags");

    NjettinessAK8Puppi_label_ = iConfig.getParameter<std::string>("NjettinessAK8Puppi_label");
    ak8PFJetsPuppi_label_     = iConfig.getParameter<std::string>("ak8PFJetsPuppi_label");


    //Consumes statements and input tags
    JetTok_                 = consumes<std::vector<pat::Jet> >       ( iConfig.getParameter<edm::InputTag>( "jetSrc"             ) );
    PuppiJetsSrc_Tok_       = consumes<std::vector<pat::Jet> >       ( iConfig.getParameter<edm::InputTag>( "puppiJetsSrc"       ) );
    EleLVec_Tok_            = consumes<std::vector<TLorentzVector> > ( iConfig.getParameter<edm::InputTag>( "eleLVec"            ) );
    MuLVec_Tok_             = consumes<std::vector<TLorentzVector> > ( iConfig.getParameter<edm::InputTag>( "muLVec"             ) );
    TrksForIsoVetolVec_Tok_ = consumes<std::vector<TLorentzVector> > ( iConfig.getParameter<edm::InputTag>( "trksForIsoVetoLVec" ) );
    LooseIsoTrksVec_Tok_    = consumes<std::vector<TLorentzVector> > ( iConfig.getParameter<edm::InputTag>( "looseisoTrksLVec"   ) );
    
    //DeepAK8
   /* 
    auto cc = consumesCollector();
    fatjetNN_ = std::make_unique<deepntuples::FatJetNN> deepntuples::FatJetNN(iConfig, cc);
    // load json for input variable transformation
    fatjetNN_->load_json("preprocessing.json"); // use the full path or put the file in the current working directory (i.e., where you run cmsRun)
    // load DNN model and parameter files
    fatjetNN_->load_model("resnet-symbol.json", "resnet.params"); // use the full path or put the file in the current working directory (i.e., where you run cmsRun)
    */
    
     auto cc = consumesCollector();
     fatjetNN_ = std::make_unique<deepntuples::FatJetNN>(iConfig, cc);
     const auto& datapath_(iConfig.getParameter<std::string>("datapath"));
     // Load json for input variable transformation
     fatjetNN_->load_json(edm::FileInPath(datapath_+"/preprocessing.json").fullPath());
     // Load DNN model and parameter files
     fatjetNN_->load_model(edm::FileInPath(datapath_+"/resnet-symbol.json").fullPath(), edm::FileInPath(datapath_+"/resnet.params").fullPath());
    
    //Produces statements 
    
    //AK4 jet variables 
    produces<std::vector<TLorentzVector> >("jetsLVec");
    produces<std::vector<int> >("recoJetsFlavor");
    produces<std::vector<float> >("recoJetsCSVv2");
    produces<std::vector<float> >("recoJetsCharge");
    produces<std::vector<float> >("recoJetsJecUnc");
    produces<std::vector<float> >("recoJetsJecScaleRawToFull");
    produces<int>("nJets");

    produces<std::vector<float> >("DeepCSVall");
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

    produces<std::vector<float> >("qgLikelihood");
    produces<std::vector<float> >("qgPtD");
    produces<std::vector<float> >("qgAxis2");
    produces<std::vector<float> >("qgAxis1");
    produces<std::vector<float> >("qgMult");

    produces<std::vector<float> >("recoJetschargedHadronEnergyFraction");
    produces<std::vector<float> >("recoJetschargedEmEnergyFraction");
    produces<std::vector<float> >("recoJetsneutralEmEnergyFraction");
    produces<std::vector<float> >("recoJetsHFHadronEnergyFraction");
    produces<std::vector<float> >("recoJetsmuonEnergyFraction");
    produces<std::vector<float> >("recoJetsneutralEnergyFraction");
    produces<std::vector<float> >("recoJetsHFEMEnergyFraction");  
    produces<std::vector<float> >("PhotonEnergyFraction");
    produces<std::vector<float> >("ElectronEnergyFraction");

    produces<std::vector<float> >("ChargedHadronMultiplicity");
    produces<std::vector<float> >("NeutralHadronMultiplicity");
    produces<std::vector<float> >("PhotonMultiplicity");
    produces<std::vector<float> >("ElectronMultiplicity");
    produces<std::vector<float> >("MuonMultiplicity");

    produces<std::vector<float> >("JetProba");
    produces<std::vector<float> >("JetBprob");

    produces<std::vector<int> >("muMatchedJetIdx");
    produces<std::vector<int> >("eleMatchedJetIdx");
    produces<std::vector<int> >("trksForIsoVetoMatchedJetIdx");
    produces<std::vector<int> >("looseisoTrksMatchedJetIdx");


    //AK8 jet variables 
    produces<std::vector<TLorentzVector> >("puppiJetsLVec");
    produces<std::vector<float> >("puppiSoftDropMass");
    produces<std::vector<float> >("puppiTau1");
    produces<std::vector<float> >("puppiTau2");
    produces<std::vector<float> >("puppiTau3");
    produces<std::vector<std::vector<TLorentzVector> > >("puppiAK8SubjetLVec");
    produces<std::vector<std::vector<float> > >("puppiAK8SubjetMult");
    produces<std::vector<std::vector<float> > >("puppiAK8SubjetPtD");
    produces<std::vector<std::vector<float> > >("puppiAK8SubjetAxis1");
    produces<std::vector<std::vector<float> > >("puppiAK8SubjetAxis2");
    produces<std::vector<std::vector<float> > >("puppiAK8SubjetBDisc");

    //DeepAK8
    produces<std::vector<float> > ("puppiDeepAK8btop");
    produces<std::vector<float> > ("puppiDeepAK8bW");
    produces<std::vector<float> > ("puppiDeepAK8bZ");
    produces<std::vector<float> > ("puppiDeepAK8bZbb");
    produces<std::vector<float> > ("puppiDeepAK8bHbb");
    produces<std::vector<float> > ("puppiDeepAK8bH4q");
    produces<std::vector<std::vector<float>> > ("puppiDeepAK8raw"); 

}


prodJets::~prodJets() 
{
}

void prodJets::produceAK4JetVariables(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
    //Get AK4 jets
    edm::Handle<std::vector<pat::Jet> > jets; 
    iEvent.getByToken(JetTok_, jets);

    //get the JEC uncertainties
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get(jetType_, JetCorParColl);
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    std::unique_ptr<JetCorrectionUncertainty> jecUnc( new JetCorrectionUncertainty(JetCorPar) );

    //get the electron and muon collections saved to the tuples
    edm::Handle<std::vector<TLorentzVector> > eleLVec, muLVec;
    iEvent.getByToken(EleLVec_Tok_, eleLVec);
    iEvent.getByToken(MuLVec_Tok_, muLVec);

    //get the isotrack collections saved to the tuples 
    edm::Handle<std::vector<TLorentzVector> > trksForIsoVetoLVec, looseisoTrksLVec;
    iEvent.getByToken(TrksForIsoVetolVec_Tok_, trksForIsoVetoLVec);
    iEvent.getByToken(LooseIsoTrksVec_Tok_,looseisoTrksLVec);

    //Declare smart pointers and output vectors 
    //Basic jet properties 
    std::unique_ptr<std::vector<TLorentzVector> > jetsLVec(new std::vector<TLorentzVector>());
    std::unique_ptr<std::vector<int> > recoJetsFlavor(new std::vector<int>());
    std::unique_ptr<std::vector<float> > recoJetsCSVv2(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsCharge(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsJecUnc(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsJecScaleRawToFull(new std::vector<float>());

    //Jet porbabilities 
    std::unique_ptr<std::vector<float> > JetProba(new std::vector<float>());
    std::unique_ptr<std::vector<float> > JetBprob(new std::vector<float>());

    //DeepCSV jet tagger
    std::unique_ptr<std::vector<float> > DeepCSVall(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepCSVb(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepCSVc(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepCSVl(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepCSVbb(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepCSVcc(new std::vector<float>());

    //DeepFlavour jet tagger
    std::unique_ptr<std::vector<float> > DeepFlavorb(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepFlavorbb(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepFlavorlepb(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepFlavorc(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepFlavoruds(new std::vector<float>());
    std::unique_ptr<std::vector<float> > DeepFlavorg(new std::vector<float>());

    //Legacy charm tagger
    std::unique_ptr<std::vector<float> > CversusB(new std::vector<float>());
    std::unique_ptr<std::vector<float> > CversusL(new std::vector<float>());

    //quark-gluon liklihood jet shape variables
    std::unique_ptr<std::vector<float> > qgLikelihood(new std::vector<float>());
    std::unique_ptr<std::vector<float> > qgPtD(new std::vector<float>());
    std::unique_ptr<std::vector<float> > qgAxis1(new std::vector<float>());
    std::unique_ptr<std::vector<float> > qgAxis2(new std::vector<float>());
    std::unique_ptr<std::vector<float> > qgMult(new std::vector<float>());

    //PF energy fraction vectors 
    std::unique_ptr<std::vector<float> > recoJetschargedHadronEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetschargedEmEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsneutralEmEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsmuonEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsneutralEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsHFEMEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > recoJetsHFHadronEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > PhotonEnergyFraction(new std::vector<float>());
    std::unique_ptr<std::vector<float> > ElectronEnergyFraction(new std::vector<float>());

    //PF object multiplicity vectors 
    std::unique_ptr<std::vector<float> > ChargedHadronMultiplicity(new std::vector<float>());
    std::unique_ptr<std::vector<float> > NeutralHadronMultiplicity(new std::vector<float>());
    std::unique_ptr<std::vector<float> > PhotonMultiplicity(new std::vector<float>());
    std::unique_ptr<std::vector<float> > ElectronMultiplicity(new std::vector<float>());
    std::unique_ptr<std::vector<float> > MuonMultiplicity(new std::vector<float>());

    //Matching vectors, for dR matching between jets and other objects  
    std::unique_ptr<std::vector<int> > muMatchedJetIdx(new std::vector<int>(muLVec->size(), -1));
    std::unique_ptr<std::vector<int> > eleMatchedJetIdx(new std::vector<int>(eleLVec->size(), -1));
    std::unique_ptr<std::vector<int> > trksForIsoVetoMatchedJetIdx(new std::vector<int>(trksForIsoVetoLVec->size(), -1));
    std::unique_ptr<std::vector<int> > looseisoTrksMatchedJetIdx(new std::vector<int>(looseisoTrksLVec->size(), -1));

    for(unsigned int ij=0; ij < jets->size(); ij++)
    {
        const pat::Jet& jet = (*jets)[ij];

        if(!(jet.pt() > AK4minPTcut_)) continue;

        //Create TLorentz vector for the jet 
        TLorentzVector perJetLVec;
        perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
        jetsLVec->push_back(perJetLVec);

        //Get DeepCSV values
        DeepCSVall  ->push_back( jet.bDiscriminator(deepCSVBJetTags_+":proball")    );
        DeepCSVb  ->push_back( jet.bDiscriminator(deepCSVBJetTags_+":probb")    );
        DeepCSVc  ->push_back( jet.bDiscriminator(deepCSVBJetTags_+":probc")    );
        DeepCSVl  ->push_back( jet.bDiscriminator(deepCSVBJetTags_+":probudsg") );
        DeepCSVbb ->push_back( jet.bDiscriminator(deepCSVBJetTags_+":probbb")   );
        DeepCSVcc ->push_back( jet.bDiscriminator(deepCSVBJetTags_+":probcc")   );

        //JetProb values
        JetProba ->push_back( jet.bDiscriminator(jetPBJetTags_)  );
        JetBprob ->push_back( jet.bDiscriminator(jetBPBJetTags_) );

        //Legacy charm tagger (for SUS-16-049 BDT resolved tagger)
        CversusB ->push_back( jet.bDiscriminator(CvsBCJetTags_) );
        CversusL ->push_back( jet.bDiscriminator(CvsLCJetTags_) );

        //deepFlavor
        DeepFlavorb    ->push_back( jet.bDiscriminator(deepFlavorBJetTags_+":probb")    );
        DeepFlavorbb   ->push_back( jet.bDiscriminator(deepFlavorBJetTags_+":probbb")   );
        DeepFlavorlepb ->push_back( jet.bDiscriminator(deepFlavorBJetTags_+":problepb") );
        DeepFlavorc    ->push_back( jet.bDiscriminator(deepFlavorBJetTags_+":probc")    );
        DeepFlavoruds  ->push_back( jet.bDiscriminator(deepFlavorBJetTags_+":probuds")  );
        DeepFlavorg    ->push_back( jet.bDiscriminator(deepFlavorBJetTags_+":probg")    );

        //Additional jec qualities
        float scaleRawToFull = jet.jecFactor(jet.currentJECLevel(), "none", jet.currentJECSet())/jet.jecFactor("Uncorrected", "none", jet.currentJECSet());
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

        //jet parton flavor 
        recoJetsFlavor ->push_back( jet.partonFlavour() );

        //quark-gluon liklihood jet shape variables 
        qgLikelihood   ->push_back( jet.userFloat(qgTaggerKey_+":qgLikelihood") );
        qgPtD          ->push_back( jet.userFloat(qgTaggerKey_+":ptD")          );
        qgAxis1        ->push_back( jet.userFloat(qgTaggerKey_+":axis1")        );
        qgAxis2        ->push_back( jet.userFloat(qgTaggerKey_+":axis2")        );
        qgMult         ->push_back( jet.userInt(  qgTaggerKey_+":mult")         );

        //CSVv2 b-tag discriminator 
        recoJetsCSVv2 ->push_back( jet.bDiscriminator(bTagKeyString_) );

        //Jet charge 
        recoJetsCharge ->push_back( jet.jetCharge() );

        //PF energy fractions for the jet 
        recoJetschargedHadronEnergyFraction ->push_back( jet.chargedHadronEnergyFraction() );
        recoJetsneutralEnergyFraction       ->push_back( jet.neutralHadronEnergyFraction() );
        recoJetschargedEmEnergyFraction     ->push_back( jet.chargedEmEnergyFraction()     );
        recoJetsneutralEmEnergyFraction     ->push_back( jet.neutralEmEnergyFraction()     );
        recoJetsmuonEnergyFraction          ->push_back( jet.muonEnergyFraction()          );
        PhotonEnergyFraction                ->push_back( jet.photonEnergyFraction()        );
        ElectronEnergyFraction              ->push_back( jet.electronEnergyFraction()      );
        recoJetsHFHadronEnergyFraction      ->push_back( jet.HFHadronEnergyFraction()      );
        recoJetsHFEMEnergyFraction          ->push_back( jet.HFEMEnergyFraction()          );

        //PF particle type multiplicities 
        ChargedHadronMultiplicity ->push_back( jet.chargedHadronMultiplicity() );
        NeutralHadronMultiplicity ->push_back( jet.neutralHadronMultiplicity() );
        PhotonMultiplicity        ->push_back( jet.photonMultiplicity()        );
        ElectronMultiplicity      ->push_back( jet.electronMultiplicity()      );
        MuonMultiplicity          ->push_back( jet.muonMultiplicity()          );
        
        //calculate matching vectors
        drMatchObject(muLVec,             jet, ij, muMatchedJetIdx);
        drMatchObject(eleLVec,            jet, ij, eleMatchedJetIdx);
        drMatchObject(trksForIsoVetoLVec, jet, ij, trksForIsoVetoMatchedJetIdx);
        drMatchObject(looseisoTrksLVec,   jet, ij, looseisoTrksMatchedJetIdx);
    }

    std::unique_ptr<int> nJets (new int(jetsLVec->size()));

    // store in the event
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
 
    iEvent.put(std::move(DeepCSVall),"DeepCSVall"); 
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

}

void prodJets::produceAK8JetVariables(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
    //PUPPI AK8 jets and subjets
    edm::Handle<std::vector<pat::Jet> > puppiJets;
    iEvent.getByToken(PuppiJetsSrc_Tok_, puppiJets);

    //deepntuples::FatJetNN *fatjetNN_;

    //Create smart pointer for output vectors 
    std::unique_ptr<std::vector<TLorentzVector> > puppiJetsLVec(new std::vector<TLorentzVector>());
    std::unique_ptr<std::vector<float> > puppiSoftDropMass(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiTau1(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiTau2(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiTau3(new std::vector<float>());

    std::unique_ptr<std::vector<std::vector<TLorentzVector> > > puppiAK8SubjetLVec(new std::vector<std::vector<TLorentzVector> >());
    std::unique_ptr<std::vector<std::vector<float> > > puppiAK8SubjetMult(new std::vector<std::vector<float > >());
    std::unique_ptr<std::vector<std::vector<float> > > puppiAK8SubjetPtD(new std::vector<std::vector<float > >());
    std::unique_ptr<std::vector<std::vector<float> > > puppiAK8SubjetAxis1(new std::vector<std::vector<float > >());
    std::unique_ptr<std::vector<std::vector<float> > > puppiAK8SubjetAxis2(new std::vector<std::vector<float > >());
    std::unique_ptr<std::vector<std::vector<float> > > puppiAK8SubjetBDisc(new std::vector<std::vector<float > >());

    //DeepAK8
    std::unique_ptr<std::vector<float> > puppiDeepAK8btop(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiDeepAK8bW(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiDeepAK8bZ(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiDeepAK8bZbb(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiDeepAK8bHbb(new std::vector<float>());
    std::unique_ptr<std::vector<float> > puppiDeepAK8bH4q(new std::vector<float>());
    std::unique_ptr<std::vector<std::vector<float>> > puppiDeepAK8raw(new std::vector<std::vector<float>>());

    fatjetNN_->readEvent(iEvent, iSetup);

    for(unsigned int ip = 0; ip < puppiJets->size(); ip++)
    {
        auto& fatjet = (*puppiJets)[ip];
        if((fatjet.pt())>AK8minPTcut_)
        {
            puppiJetsLVec->emplace_back(fatjet.p4().X(), fatjet.p4().Y(), fatjet.p4().Z(), fatjet.p4().T());

            puppiTau1 ->push_back( fatjet.userFloat( NjettinessAK8Puppi_label_+":tau1" ) );
            puppiTau2 ->push_back( fatjet.userFloat( NjettinessAK8Puppi_label_+":tau2" ) );
            puppiTau3 ->push_back( fatjet.userFloat( NjettinessAK8Puppi_label_+":tau3" ) );

            puppiSoftDropMass ->push_back( fatjet.userFloat( ak8PFJetsPuppi_label_+"SoftDropMass" ) );

            puppiAK8SubjetLVec->push_back(std::vector<TLorentzVector>());
            puppiAK8SubjetMult->push_back(std::vector<float>());
            puppiAK8SubjetPtD->push_back(std::vector<float>());
            puppiAK8SubjetAxis1->push_back(std::vector<float>());
            puppiAK8SubjetAxis2->push_back(std::vector<float>());
            puppiAK8SubjetBDisc->push_back(std::vector<float>());     

            // run the NN predictions
            deepntuples::JetHelper jet_helper(&fatjet);
            const auto& nnpreds = fatjetNN_->predict(jet_helper);
            deepntuples::FatJetNNHelper nn(nnpreds);
            
            // get the scores
            puppiDeepAK8btop->push_back(nn.get_binarized_score_top());
            puppiDeepAK8bW->push_back(nn.get_binarized_score_w());
            puppiDeepAK8bZ->push_back(nn.get_binarized_score_z());
            puppiDeepAK8bZbb->push_back(nn.get_binarized_score_zbb());
            puppiDeepAK8bHbb->push_back(nn.get_binarized_score_hbb());
            puppiDeepAK8bH4q->push_back(nn.get_binarized_score_h4q());
            puppiDeepAK8raw->push_back(std::vector<float>());
            //std::cout<< "AK8jet_PT " << fatjet.pt() << " AK8Jet eta "<< fatjet.eta() << " get binarized score "<< nn.get_binarized_score_top()<< " W score "<<nn.get_binarized_score_w() << " Zbb score " << nn.get_binarized_score_zbb() << std::endl;
            for(const auto& pred : nnpreds) puppiDeepAK8raw->back().push_back(pred);

            auto const & subjets = fatjet.subjets("SoftDrop");
            for( auto const & it : subjets)
            {
                TLorentzVector perSubJetLVec;
                perSubJetLVec.SetPtEtaPhiE( it->pt(), it->eta(), it->phi(), it->energy() );
          
                // btag info
                double subjetBDiscriminator = it->bDiscriminator(bTagKeyString_.c_str());
      
                //compute the qg input variables for the subjet
                double totalMult = 0;
                double ptD       = 0;
                double axis1     = 0;
                double axis2     = 0;
                computeJetShape(dynamic_cast<const pat::Jet *>(&(*it)), true, totalMult, ptD, axis1, axis2);
      
                puppiAK8SubjetLVec->back().push_back(perSubJetLVec);
                puppiAK8SubjetMult->back().push_back(totalMult);
                puppiAK8SubjetPtD->back().push_back(ptD);
                puppiAK8SubjetAxis1->back().push_back(axis1);
                puppiAK8SubjetAxis2->back().push_back(axis2);
                puppiAK8SubjetBDisc->back().push_back(subjetBDiscriminator);
            }
        }
    }
    //put variables back into event
    iEvent.put(std::move(puppiJetsLVec), "puppiJetsLVec");
    iEvent.put(std::move(puppiSoftDropMass), "puppiSoftDropMass");
    iEvent.put(std::move(puppiTau1),"puppiTau1");
    iEvent.put(std::move(puppiTau2),"puppiTau2");
    iEvent.put(std::move(puppiTau3),"puppiTau3");

    iEvent.put(std::move(puppiAK8SubjetLVec),  "puppiAK8SubjetLVec");
    iEvent.put(std::move(puppiAK8SubjetMult),  "puppiAK8SubjetMult");
    iEvent.put(std::move(puppiAK8SubjetPtD),   "puppiAK8SubjetPtD");
    iEvent.put(std::move(puppiAK8SubjetAxis1), "puppiAK8SubjetAxis1");
    iEvent.put(std::move(puppiAK8SubjetAxis2), "puppiAK8SubjetAxis2");
    iEvent.put(std::move(puppiAK8SubjetBDisc), "puppiAK8SubjetBDisc");

    //DeepAK8
    iEvent.put(std::move(puppiDeepAK8btop), "puppiDeepAK8btop");
    iEvent.put(std::move(puppiDeepAK8bW), "puppiDeepAK8bW");
    iEvent.put(std::move(puppiDeepAK8bZ), "puppiDeepAK8bZ");
    iEvent.put(std::move(puppiDeepAK8bZbb), "puppiDeepAK8bZbb");
    iEvent.put(std::move(puppiDeepAK8bHbb), "puppiDeepAK8bHbb");
    iEvent.put(std::move(puppiDeepAK8bH4q), "puppiDeepAK8bH4q");
    iEvent.put(std::move(puppiDeepAK8raw), "puppiDeepAK8raw");
}

bool prodJets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
    produceAK4JetVariables(iEvent, iSetup);
    produceAK8JetVariables(iEvent, iSetup);
        
    return true;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(prodJets);
