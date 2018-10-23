// -*- C++ -*-
//
// Package:    SuSySubstructure
// Class:      PhotonIDisoProducer
// 
/*

  Description: Takes as cfg input a photon collection
  recomputes sigmaIetaIeta, applies loose EGamma WP cuts,
  fills 4-vector information for the best photon, ID & ISO
  variables for all photons, and counts the number of good
  photons.
  
*/
//
// Original Author:  Andrew Whitbeck
//         Created:  Wed March 7, 2014
// 

// system include files
#include <memory>
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhotonIDisoProducer.h"
#include "effArea.cc"

#include "TLorentzVector.h"
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include <vector>

PhotonIDisoProducer::PhotonIDisoProducer(const edm::ParameterSet& iConfig):
  photonCollection(iConfig.getUntrackedParameter<edm::InputTag>("photonCollection")),
  electronCollection(iConfig.getUntrackedParameter<edm::InputTag>("electronCollection")),
  conversionCollection(iConfig.getUntrackedParameter<edm::InputTag>("conversionCollection")),
  beamspotCollection(iConfig.getUntrackedParameter<edm::InputTag>("beamspotCollection")),
  ecalRecHitsInputTag_EE_(iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE")),
  ecalRecHitsInputTag_EB_(iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB")),
  rhoCollection(iConfig.getUntrackedParameter<edm::InputTag>("rhoCollection")),
  genParCollection(iConfig.getUntrackedParameter<edm::InputTag>("genParCollection")),
  loosePhotonID(iConfig.getParameter<edm::InputTag>("loosePhotonID")),
  mediumPhotonID(iConfig.getParameter<edm::InputTag>("mediumPhotonID")),
  tightPhotonID(iConfig.getParameter<edm::InputTag>("tightPhotonID")),
  photonTok_(consumes<edm::View<pat::Photon>>(photonCollection)),
  electronTok_(consumes<pat::ElectronCollection>(electronCollection)),
  conversionTok_(consumes<std::vector<reco::Conversion>>(conversionCollection)),
  beamspotTok_(consumes<reco::BeamSpot>(beamspotCollection)),
  ecalRecHitsInputTag_EE_Token_(consumes<EcalRecHitCollection>(ecalRecHitsInputTag_EE_)),
  ecalRecHitsInputTag_EB_Token_(consumes<EcalRecHitCollection>(ecalRecHitsInputTag_EB_)),
  rhoTok_(consumes<double>(rhoCollection)),
  genParTok_(consumes<edm::View<reco::GenParticle>>(genParCollection)),
  looseIdToken_(consumes<edm::ValueMap<bool> >(loosePhotonID)),
  mediumIdToken_(consumes<edm::ValueMap<bool> >(mediumPhotonID)),
  tightIdToken_(consumes<edm::ValueMap<bool> >(tightPhotonID)),
  debug(iConfig.getUntrackedParameter<bool>("debug",true))
{

  ecalRecHitsInputTag_EE_Token_ = consumes<EcalRecHitCollection>(ecalRecHitsInputTag_EE_);
  ecalRecHitsInputTag_EB_Token_ = consumes<EcalRecHitCollection>(ecalRecHitsInputTag_EB_);

  //Andres TLorentz 
  produces<std::vector<TLorentzVector> >("photonLVec");
  produces<std::vector<TLorentzVector> >("photonLVecGen");
  produces<std::vector<TLorentzVector> >("genPartonLVec");

  produces< std::vector< pat::Photon > >(); 
  produces< std::vector< float > >("photonISEB");
  produces< std::vector< float > >("photonGenMatched"); 
  produces< std::vector< float > >("photonHadTowOverEM"); 
  produces< std::vector< float > >("photonSigmaIetaIeta"); 
  produces< std::vector< float > >("photonPFChargedIso"); 
  produces< std::vector< float > >("photonPFNeutralIso"); 
  produces< std::vector< float > >("photonPFGammaIso"); 
  produces< std::vector< float > >("photonPFChargedIsoRhoCorr"); 
  produces< std::vector< float > >("photonPFNeutralIsoRhoCorr"); 
  produces< std::vector< float > >("photonPFGammaIsoRhoCorr"); 
  produces< std::vector< float > >("photonHasPixelSeed"); 
  produces< std::vector< float > >("passElectronVeto"); 
  //produces< std::vector< bool > >("hadronization");
  produces< std::vector< bool > >("photonNonPrompt");
  produces< std::vector< bool > >("fullID");
  produces< std::vector< bool > >("loosePhotonID");
  produces< std::vector< bool > >("mediumPhotonID");
  produces< std::vector< bool > >("tightPhotonID");
  produces< std::vector< bool > >("extraLooseID");  
  //Andres: adding pt, eta and phi
  produces< std::vector< float > >("photonPt");
  produces< std::vector< float > >("photonEta");
  produces< std::vector< float > >("photonPhi");

}


PhotonIDisoProducer::~PhotonIDisoProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}
//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonIDisoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
 
  auto goodPhotons  = std::make_unique<std::vector<pat::Photon>>();
  auto photon_isEB = std::make_unique<std::vector<float>>();
  auto photon_genMatched = std::make_unique<std::vector<float>>();
  auto photon_hadTowOverEM = std::make_unique<std::vector<float>>();
  auto photon_sigmaIetaIeta = std::make_unique<std::vector<float>>();
  auto photonpfGammaIso = std::make_unique<std::vector<float>>();
  auto photon_pfChargedIso = std::make_unique<std::vector<float>>();
  auto photon_pfNeutralIso = std::make_unique<std::vector<float>>();
  auto photonpfGammaIsoRhoCorr = std::make_unique<std::vector<float>>();
  auto photon_pfChargedIsoRhoCorr = std::make_unique<std::vector<float>>();
  auto photon_pfNeutralIsoRhoCorr = std::make_unique<std::vector<float>>();
  auto photon_hasPixelSeed = std::make_unique<std::vector<float>>();
  auto photon_passElectronVeto = std::make_unique<std::vector<float>>();
  auto   photon_hadronization = std::make_unique<std::vector<bool>>();
  auto   photon_nonPrompt  = std::make_unique<std::vector<bool>>();
  auto   photon_fullID  = std::make_unique<std::vector<bool>>();
  auto   photon_looseID  = std::make_unique<std::vector<bool>>();
  auto   photon_mediumID  = std::make_unique<std::vector<bool>>();
  auto   photon_tightID  = std::make_unique<std::vector<bool>>();
  auto   photon_electronFakes  = std::make_unique<std::vector<bool>>(); 
  auto photon_extraLooseID = std::make_unique<std::vector<bool>>();
  auto photonLVec = std::make_unique<std::vector<TLorentzVector> >();
  auto photonLVecGen = std::make_unique<std::vector<TLorentzVector> >();
  auto genPartonLVec = std::make_unique<std::vector<TLorentzVector> >();
  auto photon_pt = std::make_unique<std::vector<float>>();
  auto photon_eta = std::make_unique<std::vector<float>>();
  auto photon_phi = std::make_unique<std::vector<float>>();
  

  if( debug ){
    std::cout << "new events" << std::endl;
    std::cout << "===================" << std::endl;
  }

  Handle< View< pat::Photon> > photonCands;
  iEvent.getByToken( photonTok_,photonCands);
  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronTok_, electrons);
  Handle<vector<reco::Conversion> > conversions;
  iEvent.getByToken(conversionTok_,conversions);
  Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamspotTok_,beamSpot);
  edm::Handle< View<reco::GenParticle> > genParticles;
  iEvent.getByToken( genParTok_,genParticles);
  edm::Handle< double > rho_;
  iEvent.getByToken(rhoTok_,rho_);
  double rho = *rho_;

  edm::Handle<edm::ValueMap<bool> >   loose_id_decisions_;
  iEvent.getByToken(looseIdToken_,loose_id_decisions_);
  edm::Handle<edm::ValueMap<bool> >   medium_id_decisions_;
  iEvent.getByToken(mediumIdToken_,medium_id_decisions_);
  edm::Handle<edm::ValueMap<bool> >   tight_id_decisions_;
  iEvent.getByToken(tightIdToken_,tight_id_decisions_);

  
  //Gen level Photons
  if (genParticles.isValid()){
    for(unsigned int ig = 0; ig < genParticles->size(); ig++) {

      TLorentzVector genPhoton;
      TLorentzVector genParton;

      if (genParticles->at(ig).pt() > 10) {
	if( genParticles->at(ig).pdgId() == 22 && (genParticles->at(ig).status() == 1 || genParticles->at(ig).status() == 2 
						   || (genParticles->at(ig).status())/10 == 2)){

	  genPhoton.SetPtEtaPhiE( genParticles->at(ig).pt(), genParticles->at(ig).eta(), genParticles->at(ig).phi(), genParticles->at(ig).energy() );

	  photonLVecGen->push_back(genPhoton);

	  }
        if ((abs(genParticles->at(ig).pdgId()) == 1 || abs(genParticles->at(ig).pdgId()) == 2 || abs(genParticles->at(ig).pdgId()) == 3 ||
             abs(genParticles->at(ig).pdgId()) == 4 || abs(genParticles->at(ig).pdgId()) == 5 || abs(genParticles->at(ig).pdgId()) == 6 ||
             abs(genParticles->at(ig).pdgId()) == 9 || abs(genParticles->at(ig).pdgId()) == 21) && (genParticles->at(ig).status() == 23 ||
                                                                                          genParticles->at(ig).status() == 71)){

          genParton.SetPtEtaPhiE( genParticles->at(ig).pt(), genParticles->at(ig).eta(), genParticles->at(ig).phi(), genParticles->at(ig).energy() );
          genPartonLVec->push_back(genParton);
      }
    }
  }
 }//ADDed


  // - - - - - - - - - - - - - - - - - - - - 
  // Initializing effective area to be used 
  // for rho corrections to the photon isolation
  // variables. 
  // - - - - - - - - - - - - - - - - - - - - 
  effArea effAreas;
  effAreas.addEffA( 0.0, 1.0, 0.0234, 0.0053, 0.078 );
  effAreas.addEffA( 1.0, 1.479, 0.0189, 0.0130, 0.0629 );
  effAreas.addEffA( 1.479, 2.0, 0.0171, 0.0057, 0.0264 );
  effAreas.addEffA( 2.0, 2.2, 0.0129, 0.0070, 0.0462 );
  effAreas.addEffA( 2.2, 2.3, 0.0110, 0.0152, 0.0740 );
  effAreas.addEffA( 2.3, 2.4, 0.0074, 0.0232, 0.0924 );
  effAreas.addEffA( 2.4, 99., 0.0035, 0.1709, 0.1484 );

  /// setup cluster tools
  noZS::EcalClusterLazyTools clusterTools_(iEvent, iSetup, ecalRecHitsInputTag_EB_Token_, ecalRecHitsInputTag_EE_Token_);
        
  for( View< pat::Photon >::const_iterator iPhoton = photonCands->begin(); iPhoton != photonCands->end(); ++iPhoton){

    if( debug ) {
      std::cout << "photon pt: " << iPhoton->pt() << std::endl; 
      std::cout << "photon eta: " << iPhoton->eta() << std::endl;
      std::cout << "photon phi: " << iPhoton->phi() << std::endl;
    }
    //ID decisions
    const edm::Ptr<pat::Photon> phoPtr(photonCands, iPhoton - photonCands->begin() );
    bool passloose = (*loose_id_decisions_)[ phoPtr ];
    bool passmedium = (*medium_id_decisions_)[ phoPtr ];
    bool passtight = (*tight_id_decisions_)[ phoPtr ];

    std::vector<float> vCov = clusterTools_.localCovariances( *(iPhoton->superCluster()->seed()) ); 
    const float sieie = (isnan(vCov[0]) ? 0. : sqrt(vCov[0])); 
    
    double chIso = effAreas.rhoCorrectedIso(  pfCh  , iPhoton->chargedHadronIso() , iPhoton->eta() , rho ); 
    double nuIso = effAreas.rhoCorrectedIso(  pfNu  , iPhoton->neutralHadronIso() , iPhoton->eta() , rho ); 
    double gamIso = effAreas.rhoCorrectedIso( pfGam , iPhoton->photonIso()        , iPhoton->eta() , rho ); 

    // apply photon selection -- all good photons will be saved
    // use loose selection with no sieie or chiso cuts
    bool isBarrelPhoton=false;
    bool isEndcapPhoton=false;
    bool passAcc=false;

    double PhEta=iPhoton->eta();

    if(fabs(PhEta) < 1.4442  ){
         isBarrelPhoton=true;
                    }
    
    else if(fabs(PhEta)>1.566 && fabs(PhEta)<2.5){
                    isEndcapPhoton=true;
           }
    else {
        isBarrelPhoton=false;
        isEndcapPhoton=false;
        }
    
     if(isBarrelPhoton || isEndcapPhoton){
              passAcc=true;
            }

    // check if photon is a good loose photon
    if( passAcc ){ 
      TLorentzVector perGammaLVec; 
      perGammaLVec.SetPtEtaPhiE( iPhoton->pt(), iPhoton->eta(), iPhoton->phi(), iPhoton->energy() );
      photonLVec->push_back(perGammaLVec);
      
      goodPhotons->push_back( *iPhoton );
      photon_isEB->push_back( iPhoton->isEB() );
      photon_genMatched->push_back( iPhoton->genPhoton() != NULL );
      photon_hadTowOverEM->push_back( iPhoton->hadTowOverEm() ) ;
      photon_sigmaIetaIeta->push_back( sieie );
      photon_pfChargedIso->push_back(      iPhoton->chargedHadronIso() );
      photonpfGammaIso->push_back(        iPhoton->photonIso() );
      photon_pfNeutralIso->push_back(      iPhoton->neutralHadronIso() );
      photon_pfChargedIsoRhoCorr->push_back( chIso  );
      photonpfGammaIsoRhoCorr->push_back(   gamIso  );
      photon_pfNeutralIsoRhoCorr->push_back( nuIso );
      photon_hasPixelSeed->push_back( iPhoton->hasPixelSeed() );
      photon_passElectronVeto->push_back( !hasMatchedPromptElectron(iPhoton->superCluster(),electrons, conversions, beamSpot->position()) );
      photon_pt->push_back( iPhoton->pt() );
      photon_eta->push_back( iPhoton->eta() );
      photon_phi->push_back( iPhoton->phi() );
      photon_looseID->push_back(passloose);
      photon_mediumID->push_back(passmedium);
      photon_tightID->push_back(passtight);
      }
      
      if (genParticles.isValid()){//genLevel Stuff
	// loop over gen particles and find nonprompt and hadronization photons
	int matchedGenPrompt = 0;
	int matchedGenNonPrompt = 0 ;
        
	for( View<reco::GenParticle>::const_iterator iGen = genParticles->begin();
	     iGen != genParticles->end();
	     ++iGen){
	  
	  // check for non-prompt photons ----------------------
	  if( iGen->pdgId() == 22 && ( ( iGen->status() / 10 ) == 2 || iGen->status() == 1 || iGen->status() == 2 ) ){
        
	    TLorentzVector gen( iGen->px() , iGen->py() , iGen->pz() , iGen->energy() );
	    TLorentzVector photon( iPhoton->px() , iPhoton->py() , iPhoton->pz() , iPhoton->energy() );
	  
	    if( gen.DeltaR(photon) < 0.2 ){ /// I LEFT OFF HERE!!!!!!
	      if( abs(iGen->mother()->pdgId()) > 100 && abs(iGen->mother()->pdgId()) != 2212 ) matchedGenNonPrompt++ ;
	      if( abs(iGen->mother()->pdgId()) <= 22 || abs(iGen->mother()->pdgId()) == 2212 ) matchedGenPrompt++ ;
            }
          }
	  
          // ----------------------------------------------------
	  
	  
        }// end loop over gen particles

	if( matchedGenPrompt > 0 || matchedGenNonPrompt == 0 ) photon_nonPrompt->push_back(false);
	else if( matchedGenNonPrompt > 0 ) photon_nonPrompt->push_back(true);
	else photon_nonPrompt->push_back(false);
      }//gen level stuff
      //photon_hadronization->push_back( isHadronization );
      
      //}//pure photons
      
  }// end loop over candidate photons
  
  iEvent.put(std::move(photon_isEB) , "photonISEB" );
  iEvent.put(std::move(photon_genMatched) , "photonGenMatched" );
  iEvent.put(std::move(photon_hadTowOverEM) , "photonHadTowOverEM" );
  iEvent.put(std::move(photon_sigmaIetaIeta) , "photonSigmaIetaIeta" );
  iEvent.put(std::move(photon_pfChargedIso) , "photonPFChargedIso" );
  iEvent.put(std::move(photon_pfNeutralIso) , "photonPFNeutralIso" );
  iEvent.put(std::move(photonpfGammaIso) , "photonPFGammaIso" );
  iEvent.put(std::move(photon_pfChargedIsoRhoCorr) , "photonPFChargedIsoRhoCorr" );
  iEvent.put(std::move(photon_pfNeutralIsoRhoCorr) , "photonPFNeutralIsoRhoCorr" );
  iEvent.put(std::move(photonpfGammaIsoRhoCorr) , "photonPFGammaIsoRhoCorr" );
  iEvent.put(std::move(photon_hasPixelSeed) , "photonHasPixelSeed" );
  iEvent.put(std::move(photon_nonPrompt) , "photonNonPrompt" );
  iEvent.put(std::move(photon_looseID) , "loosePhotonID" );
  iEvent.put(std::move(photon_mediumID) , "mediumPhotonID" );
  iEvent.put(std::move(photon_tightID) , "tightPhotonID" );
  
  //Andres
  //Gamma TLorentz
  iEvent.put(std::move(photonLVec),"photonLVec");
  iEvent.put(std::move(photonLVecGen),"photonLVecGen");
  iEvent.put(std::move(genPartonLVec), "genPartonLVec");
}

// copied from https://github.com/RazorCMS/SUSYBSMAnalysis-RazorTuplizer/blob/6072ffb43bbeb3f6b34cf8a96426c7f104c5b902/plugins/RazorAux.cc#L127
//check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
//and not matching any conversion in the collection passing the quality cuts
bool PhotonIDisoProducer::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<pat::Electron> > &eleCol,
                                                   const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,
                                                   float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax)
{
  if (sc.isNull()) return false;
  for (std::vector<pat::Electron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
    //match electron to supercluster
    if (it->superCluster()!=sc) continue;
    //check expected inner hits 2017/2018 numberofAllHits named CHANGED!!!!
    if (it->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
    //check if electron is matching to a conversion
    if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
    return true;
  }
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 

PhotonIDisoProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonIDisoProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
PhotonIDisoProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PhotonIDisoProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PhotonIDisoProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PhotonIDisoProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonIDisoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {


}


#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonIDisoProducer);
