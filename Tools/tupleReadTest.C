#include "NTupleReader.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

namespace joeFunctions
{
    void calcJoe(NTupleReader& tr)
    {
        double retval = tr.getVar<double>("met") + tr.getVar<double>("mht");

        int five = 5;

        std::vector<double> *vec = new std::vector<double>(4, tr.getVar<double>("mht"));

        tr.registerDerivedVar("joe", retval);
        tr.registerDerivedVar("five", five);
        tr.registerDerivedVec("threeNum", vec);
    }
}

int main()
{
    //char nBase[] = "root://cmsxrootd-site.fnal.gov//store/user/pastika/DYJetsToLL_M-50_13TeV-madgraph-pythia8/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141227_223539/0000/stopFlatNtuples_%d.root";
    //char nBase[] = "/eos/uscms/store/user/pastika/DYJetsToLL_M-50_13TeV-madgraph-pythia8/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141227_223539/0000/stopFlatNtuples_%d.root";
    char nBase[] = "/eos/uscms/store/user/lpcsusyhad/PHYS14_720_Mar14_2014_v2/pastika/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/150328_003540/0000/stopFlatNtuples_%d.root";
    //char nBase[] = "root://cmsxrootd-site.fnal.gov//store/user/lpcsusyhad/PHYS14_720_Dec23_2014/pastika/DYJetsToLL_M-50_13TeV-madgraph-pythia8/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141227_223539/0000/stopFlatNtuples_%d.root";
    char fb[] = "../SkimsAUX/workdir/prod/stopFlatNtuples.root";
    TChain *f = new TChain("stopTreeMaker/AUX");

    size_t t0 = clock();

    //char fname[512];
    //for(int i = 1; i <= 1; ++i)
    //{
    //    sprintf(fname, nBase, i);
    //    f->Add(fname);
    //}
    f->Add(fb);
    
    NTupleReader tr(f);
    tr.registerFunction(&joeFunctions::calcJoe);

    std::cout << "NEVTS: " << tr.getNEntries() << std::endl;

    while(tr.getNextEvent())
    {
        if(tr.getEvtNum() == 1)
        {
            tr.printTupleMembers();
            FILE * fout = fopen("NTupleTypes.txt", "w");
            tr.printTupleMembers(fout);
            fclose(fout);
        }
        if(tr.getEvtNum()%100000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
        
        const std::vector<TLorentzVector>& jetsLVec             = tr.getVec<TLorentzVector>("jetsLVec");
        const std::vector<TLorentzVector>& prodJetsJOE_jetsLVec = tr.getVec<TLorentzVector>("prodJetsJOE_jetsLVec");
        const std::vector<TLorentzVector>& muonsLVec            = tr.getVec<TLorentzVector>("muonsLVec");

        for(auto& jet : jetsLVec) std::cout << jet.Pt() << "\t";
        std::cout << std::endl;
        for(auto& jet : prodJetsJOE_jetsLVec) std::cout << jet.Pt() << "\t";
        std::cout << std::endl;
        for(auto& muon : muonsLVec) std::cout << muon.Pt() << "\t";
        std::cout << std::endl;
        std::cout << std::endl;
    }

    //const unsigned int& run = tr.getVar<unsigned int>("run");

}
