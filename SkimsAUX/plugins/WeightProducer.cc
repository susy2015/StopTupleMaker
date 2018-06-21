// -*- C++ -*-
//
// Package:    WeightProducer
// Class:      WeightProducer
//
/**\class WeightProducer WeightProducer.cc RA2/WeightProducer/src/WeightProducer.cc
 
 Description: <one line class summary>
 
 Implementation:
 <Notes on implementation>
 */
//
//         Created:  Tue Nov 10 17:58:04 CET 2009
// $Id: WeightProducer.cc,v 1.7 2012/10/31 04:16:11 lhx Exp $
//
//


// system include files
#include <cmath>
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"

//
// class decleration
//

class WeightProducer: public edm::EDProducer {
public:
    explicit WeightProducer(const edm::ParameterSet&);
    ~WeightProducer();

private:
    //  virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    const std::string _weightingMethod;

    float _expo;
// This is now used for total weight
    float _weight;
// This one is to record the weight from ptHat calculation
// We'd expect NOT to put all the information into "produces", 
// so we expect to check the configuration files for any other information
// on how the weight is produced.
    float ptHatweight_;
    float _startWeight;
    float _LumiScale;
    float _xs;
    float _NumberEvents;
    float _lumi;
    edm::InputTag _weightName;

    float weightFactor;

    int maxInitPrintOut_;
    int cntPrintOut_;
    float accuWeight_;

    float ptHat_;

    float weightWARNingLowThreshold_;
    float weightWARNingUpThreshold_;
    float ratioWARNingThreshold_;

// Default is an empty string which means there is no PU weighting needed 
// (Any of the input string is empty means no PU weighting...)
    const std::string inputPUfileMC_, inputPUfileData_;
    const std::string inputPUhistMC_, inputPUhistData_;
// This key string is the function name which does the PU re-weighting, i.e., weight3BX, weight3D, weight...
// TODO: Default is weight3D for now!
    const std::string reWeightPUmethod_;
    bool doPUweighting;

    edm::LumiReWeighting LumiWeights_;
    edm::Lumi3DReWeighting Lumi3DWeights_;
    float initPUscaleFactor_;
    bool initPUinput();
    float puweight_;
    edm::InputTag genSrc_;
    edm::EDGetTokenT<GenEventInfoProduct> GenTok_;
    edm::EDGetTokenT<float> WeightNameTok_;
};

WeightProducer::WeightProducer(const edm::ParameterSet& iConfig):
        _weightingMethod(iConfig.getParameter<std::string> ("Method")),
        _expo(iConfig.getParameter<double> ("Exponent")),
        _startWeight(iConfig.getParameter<double> ("weight")),
        _LumiScale(iConfig.getParameter<double> ("LumiScale")),
        _xs(iConfig.getParameter<double> ("XS")),
        _NumberEvents(iConfig.getParameter<double> ("NumberEvts")),
        _lumi(iConfig.getParameter<double> ("Lumi")),
        _weightName(iConfig.getParameter<edm::InputTag> ("weightName")),
        weightWARNingLowThreshold_(iConfig.getParameter<double> ("weightWARNingLowThreshold")),
        weightWARNingUpThreshold_(iConfig.getParameter<double> ("weightWARNingUpThreshold")),
        ratioWARNingThreshold_(iConfig.getParameter<double> ("ratioWARNingThreshold")),
        inputPUfileMC_(iConfig.getUntrackedParameter<std::string> ("inputPUfileMC", "")),
        inputPUfileData_(iConfig.getUntrackedParameter<std::string> ("inputPUfileData", "")),
        inputPUhistMC_(iConfig.getUntrackedParameter<std::string> ("inputPUhistMC", "pileup")),
        inputPUhistData_(iConfig.getUntrackedParameter<std::string> ("inputPUhistData", "pileup")),
        reWeightPUmethod_(iConfig.getUntrackedParameter<std::string> ("reWeightPUmethod", "weight3D")),
        initPUscaleFactor_(iConfig.getUntrackedParameter<double> ("initPUscaleFactor", 1.0)){

    if (_weightingMethod == "PtHat") {
       weightFactor = _lumi * _xs / (_NumberEvents * pow(15., _expo));
    }else{
       weightFactor = 1.0;
    }
 
     genSrc_ = iConfig.getParameter<edm::InputTag>("generatorSource");
     GenTok_ = consumes<GenEventInfoProduct>(genSrc_);    
     WeightNameTok_ = consumes<float>(_weightName);
    doPUweighting = initPUinput();

    maxInitPrintOut_ = 10;
    cntPrintOut_ = 0;
    accuWeight_ = 0;

    //register your products
    produces<float> ("weight");
    produces<float> ("ptHat");
    produces<float> ("ptHatweight");
    produces<float> ("puweight");
}

WeightProducer::~WeightProducer() {
   printf("\n==> WeightProducer   cntPrintOut_ : %10d  accuWeight_ : %10.8f  avgWeight_ : %10.8f\n\n", cntPrintOut_, accuWeight_, accuWeight_/cntPrintOut_);
}

// ------------ method called to produce the data  ------------
void WeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

// PU reweighting is only on MC
    if( iEvent.isRealData() ) doPUweighting = false;

    //---------------------------------------------------------------------------
    ///Put here code that computes a weight directly from iEvent/iSetup
    ///informations, i.e. using sample ID, pt-hat-information, etc;

    //Option 3: pt-hat weighting
    float storedWeight = 0.0;
    float ratioPtHat = -1.;
    ptHat_ = 0.0; ptHatweight_ = 1.0;
    if( !iEvent.isRealData() ){
        edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
        iEvent.getByToken(GenTok_, genEvtInfoHandle);
        if (genEvtInfoHandle.isValid()) {
            ptHat_ = ( genEvtInfoHandle->hasBinningValues() ? (genEvtInfoHandle->binningValues())[0] : 0.0);
            storedWeight = genEvtInfoHandle->weight();
        }
    }
    if (_weightingMethod == "PtHat") {
        // Get pthat
        edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
        iEvent.getByToken(GenTok_, genEvtInfoHandle);
        if ( !iEvent.isRealData() ) {
            _weight = weightFactor * pow(ptHat_, _expo);
            if( storedWeight != 0 ) ratioPtHat = _weight/(_lumi*_xs/_NumberEvents*storedWeight);
        } else {
            std::cout<<"WARNING:: PtHat information needed but not available: set weight to 1!"<<std::endl;
            _weight = 1;
        }
        ptHatweight_ = _weight;
    } 
    //Option 2: existing weight variable in the event named as in _weightName
    else if (_startWeight >= 0) {
        _weight = _startWeight;
    }    
    //Option 1: weight constant, as defined in cfg file
    else {
        edm::Handle<float> event_weight;
        iEvent.getByToken(WeightNameTok_, event_weight);
        _weight = (event_weight.isValid() ? (*event_weight) : 1.0);
    }

    ///This is to consider the lumi-uncertainty, i.e. to scale the weights up- or down by 1sigma of the lumi-scale
    ///uncertainty. In general the scale is 1.0!
    _weight *= _LumiScale;
    //std::cout << "_weight: " << _weight << std::endl;

    ///Also, here one could define look-up tables for all used samples.
    ///An identifying 'string' for the currently used sample could be read
    ///from the config file. Perhaps this can be obtained dirtectly from crab?
    //---------------------------------------------------------------------------

    puweight_ = 1.0;
    if( doPUweighting ){
       if( reWeightPUmethod_ == "weight" ){
          puweight_ = LumiWeights_.weight( iEvent );
       }else if( reWeightPUmethod_ == "weight3BX" ){
          std::cout<<"NO support of weight3BX anymore in 53X!"<<std::endl;
//          puweight_ = LumiWeights_.weight3BX( iEvent );
       }else if ( reWeightPUmethod_ == "weight3D" ){
          puweight_ = Lumi3DWeights_.weight3D( iEvent );
       }
    }

    _weight *= puweight_;

// Specially treatment of the negative storedWeight for MC. Default (for data) storedWeight is 0.
    if( storedWeight < 0 ) _weight *= -1;

    // put weight into the Event
    std::unique_ptr<float> weightPtr(new float(_weight));
    iEvent.put(std::move(weightPtr), "weight");

    std::unique_ptr<float> ptHatPtr(new float(ptHat_));
    iEvent.put(std::move(ptHatPtr), "ptHat");

    std::unique_ptr<float> ptHatweightPtr(new float(ptHatweight_));
    iEvent.put(std::move(ptHatweightPtr), "ptHatweight");

    std::unique_ptr<float> puweightPtr(new float(puweight_));
    iEvent.put(std::move(puweightPtr), "puweight");

    cntPrintOut_ ++; 
    accuWeight_ += _weight; 

// Specific monitoring of PtHat case
    if( _weightingMethod == "PtHat" ){
       if( ptHatweight_ <= weightWARNingLowThreshold_ || ptHatweight_ >= weightWARNingUpThreshold_ ){
          printf("\nWARNing  ==> WeightProducer   cntPrintOut_ : %8d  ptHatweight_ : %10.8e  ptHat_ : %10.8e  out-of-threshold of (low/up) : %10.8e/%10.8e  --> accuWeight_ : %10.8e  avgWeight_ : %10.8e\n", cntPrintOut_, ptHatweight_, ptHat_, weightWARNingLowThreshold_, weightWARNingUpThreshold_, accuWeight_, accuWeight_/cntPrintOut_);
       }
 
       if( fabs(ratioPtHat - 1.0) > ratioWARNingThreshold_ ){
          printf("\nWARNing  ==> WeightProducer   cntPrintOut_ : %8d  ptHatweight_ : %10.8e  storedWeight : %10.8e  ratioPtHat : %10.8f\n", cntPrintOut_, ptHatweight_, storedWeight, ratioPtHat);
       }
    }

    if( cntPrintOut_ < maxInitPrintOut_ ){
       printf("==> WeightProducer   cntPrintOut_ : %2d  _weight : %10.8f  accuWeight_ : %10.8f  avgWeight_ : %10.8f  --> ptHatweight_ : %10.8f  storedWeight : %10.8f  ratioPtHat : %10.8f  --> puweight : %10.8f\n", cntPrintOut_, _weight, accuWeight_, accuWeight_/cntPrintOut_, ptHatweight_, storedWeight, ratioPtHat, puweight_);
    }

}

bool WeightProducer::initPUinput(){

// No input PU root files, so don't do any PU reweighting
   if( inputPUfileMC_ =="" || inputPUfileData_ == "" ) return false;

   const char *pch =0;

// Only accept PU weighting in histograms in a root file
// For cases that one has access of text file, it is easy to fill a histogram and put it into a root file.
   if( !( (pch=strstr(inputPUfileMC_.c_str(), ".root")) && strlen(pch) ==5 ) ){
      std::cout<<"inputPUfileMC_ : "<<inputPUfileMC_.c_str()<<"  is NOT a root file. Only support PU distribution in root format! If not, please make one first!"<<std::endl;      
      return false;
   }
// Must have two root files, one is for MC and the other is for data
   if( !( (pch=strstr(inputPUfileData_.c_str(), ".root")) && strlen(pch) ==5 ) ){
      std::cout<<"inputPUfileData_ : "<<inputPUfileData_.c_str()<<"  is NOT a root file. Only support PU distribution in root format! If not, please make one first!"<<std::endl;      
      return false;
   }

// Initialize the LumiReWeighting or Lumi3DReWeighting class
   if( reWeightPUmethod_ == "weight3BX" || reWeightPUmethod_ == "weight" ){
      LumiWeights_ = edm::LumiReWeighting(inputPUfileMC_.c_str(), inputPUfileData_.c_str(), inputPUhistMC_.c_str(), inputPUhistData_.c_str());
   }else if( reWeightPUmethod_ == "weight3D" ){
//      Lumi3DWeights_ = edm::Lumi3DReWeighting(inputPUfileMC_.c_str(), inputPUfileData_.c_str(), inputPUhistMC_.c_str(), inputPUhistData_.c_str());
      Lumi3DWeights_ = edm::Lumi3DReWeighting();
      Lumi3DWeights_.weight3D_set(inputPUfileMC_.c_str(), inputPUfileData_.c_str(), inputPUhistMC_.c_str(), inputPUhistData_.c_str(), "weight3D.root");
      Lumi3DWeights_.weight3D_init(initPUscaleFactor_);
   }else{
      std::cout<<"reWeightPUmethod_ : "<<reWeightPUmethod_.c_str()<<"  is NOT supported!"<<std::endl;
      return false;
   }

//   if( pch ) delete pch;

   return true;
}

// ------------ method called once each job just before starting event loop  ------------
//void
//WeightProducer::beginJob()
//{
//}

// ------------ method called once each job just after ending the event loop  ------------
void WeightProducer::endJob() { }

//define this as a plug-in
DEFINE_FWK_MODULE( WeightProducer);
