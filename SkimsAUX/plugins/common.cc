#include "StopTupleMaker/SkimsAUX/plugins/common.h"

namespace commonFunctions
{

  double GetMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
				 const reco::Candidate* ptcl, std::string type, 
				 double rho, std::vector<double> EAValues, std::vector<double> EAEtaValues,
			         bool computeMT2Activity,
				 double r_iso_min, double r_iso_max, double kt_scale,
				 bool useEAcorr, bool charged_only) {
    
    if (ptcl->pt()<5.) return 99999.;

    if (EAValues.size()-EAEtaValues.size()!=1) std::cout << "Effective area vector size should be checked." << std::endl;

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(type=="electron") {
      if (isEndCapEle(*ptcl)) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(type=="muon") {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
    } else {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.); 
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(type=="electron") ptThresh = 0;
    double r_iso = std::max(r_iso_min, std::min(r_iso_max, kt_scale/ptcl->pt()));
    
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      double dr = deltaR(pfc, *ptcl);
      if (computeMT2Activity) {
	if (dr < r_iso || dr > 0.4) continue; // activity annulus
      }
      else if (dr > r_iso) continue;

      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(!computeMT2Activity && dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(!computeMT2Activity && dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(!computeMT2Activity && dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
	if (pfc.pt()>ptThresh){
          if(!computeMT2Activity && dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    double iso(0.);
    double EA(0.);
    double eta(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      if (useEAcorr) {
	// 
	// Eta
	eta = getEta(*ptcl);
	// EA
	EA = GetEA(eta,EAValues,EAEtaValues);
	double correction = rho * EA * (r_iso/0.3) * (r_iso/0.3);
	iso -= correction;
      }
      else iso -= 0.5*iso_pu; // DB correction
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }

    /*
    std::cout << type << " pt,eta: " << ptcl->pt() << " " << eta
	      << " chargediso,neutraliso,photoniso,rho,EA_miniiso,computeMT2Activity: "
	      << iso_ch << " "
	      << iso_nh << " "
	      << iso_ph << " "
	      << rho << " "
	      << EA << " "
	      << computeMT2Activity
	      << std::endl;
    */

    iso = iso/ptcl->pt();

    return iso;
  }

  
  double GetRA2Activity(edm::Handle<pat::JetCollection> jets, const reco::Candidate* ptcl, const bool useEME) {

    double activity=0;
    for (unsigned int ijet(0); ijet < jets->size(); ijet++){
       const pat::Jet &thisJet = (*jets)[ijet];
       if(deltaR(ptcl->eta(),ptcl->phi(),thisJet.eta(),thisJet.phi())>1.0) continue;
       double chFrac = thisJet.chargedHadronEnergy()/thisJet.energy();
       if (useEME) chFrac+=thisJet.chargedEmEnergy()/thisJet.energy();
       activity+=thisJet.pt() * chFrac;
    }

    return activity;
  }

  double GetTrackActivity(edm::Handle<pat::PackedCandidateCollection> other_pfcands, const pat::PackedCandidate* track) {
//    if (track->pt()<5.) return 99999.;
    if (track->pt()<5.) return -1.0;
    double trkiso(0.); 
    double r_iso = 0.3;
    for (const pat::PackedCandidate &other_pfc : *other_pfcands) {
        if (other_pfc.charge()==0) continue;
        double dr = deltaR(other_pfc, *track);
        if (dr < r_iso || dr > 0.4) continue; // activity annulus
        float dz_other = other_pfc.dz();
        if( fabs(dz_other) > 0.1 ) continue;
        trkiso += other_pfc.pt();
      }
      double activity = trkiso/track->pt();
      return activity;
  }

  double getEta(const reco::Candidate& obj) {
    try{
      const pat::Electron& electron = dynamic_cast<const pat::Electron&>(obj);
      return electron.superCluster()->eta();
    }catch (const std::bad_cast& e){
      return obj.eta();
    }
  }

  bool isEndCapEle(const reco::Candidate& obj) {
    try{
      const pat::Electron& electron = dynamic_cast<const pat::Electron&>(obj);
      return std::abs(electron.superCluster()->eta()) > 1.479;
    }catch (const std::bad_cast& e){
      return false;
    }
  }

}
