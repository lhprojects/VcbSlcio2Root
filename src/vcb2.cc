#include <vcb2.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/LCRelation.h>
#include <EVENT/Vertex.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h>
#include <sstream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"
#include <UTIL/PIDHandler.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

using namespace std;

vcb2 a_vcb2_instance;


vcb2::vcb2()
: Processor("vcb2"),
_output(0)
{
    _description = "Print MC Truth" ;
    
    _treeFileName="MCTruth.root";
    registerProcessorParameter( "TreeOutputFile" ,
                               "The name of the file to which the ROOT tree will be written" ,
                               _treeFileName ,
                               _treeFileName);
    
    _treeName="Tau";
    registerProcessorParameter( "TreeName" ,
                               "The name of the ROOT tree" ,
                               _treeName ,
                               _treeName);
    
//    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
//                             "MCPSIMUFSP",
//                             "the final state particle after simulation",
//                             _outmcpsimufsp,
//                             std::string("ReconstructedParticle") );
//    
//    
//    registerOutputCollection( LCIO::LCRELATION,
//                             "mcpsimurelation",
//                             " relation between MCP and Reco after simulation",
//                             _outMCPSIMURelation,
//                             std::string("RelationusedSed") );
    
    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}

// is B-bar hadron
int ISBbar(int mc){
    if(mc == 511 || mc == 521 || mc == 10511 || mc == 10521 || mc == 513 || mc == 523 || mc == 10513 || mc == 10523 || mc == 20513 || mc == 20523 || mc == 515 || mc == 525 || mc == 531 || mc == 10531 || mc == 533 || mc == 10533 || mc == 20533 || mc == 535 || mc == 541 || mc == 10541 || mc == 543 || mc == 10543 || mc == 20543 || mc == 545 || mc == -5122 || mc == -5112 || mc == -5212 || mc == -5222 || mc == -5114 || mc == -5214 || mc == -5224 || mc == -5132 || mc == -5232 || mc == -5312 || mc == -5322 || mc == -5314 || mc == -5324 || mc == -5332 || mc == -5334 || mc == -5142 || mc == -5242 || mc == -5412 || mc == -5422 || mc == -5414 || mc == -5424 || mc == -5342 || mc == -5432 || mc == -5434 || mc == -5442 || mc == -5444 || mc == -5512 || mc == -5522 || mc == -5514 || mc == -5524 || mc == -5532 || mc == -5534 || mc == -5542 || mc == -5544 || mc == -5544 ) {return 1;}
    else {return 0;}
}

// Is B hadron
int ISB(int mc){
    if(mc == -511 || mc == -521 || mc == -10511 || mc == -10521 || mc == -513 || mc == -523 || mc == -10513 || mc == -10523 || mc == -20513 || mc == -20523 || mc == -515 || mc == -525 || mc == -531 || mc == -10531 || mc == -533 || mc == -10533 || mc == -20533 || mc == -535 || mc == -541 || mc == -10541 || mc == -543 || mc == -10543 || mc == -20543 || mc == -545 || mc == 5122 || mc == 5112 || mc == 5212 || mc == 5222 || mc == 5114 || mc == 5214 || mc == 5224 || mc == 5132 || mc == 5232 || mc == 5312 || mc == 5322 || mc == 5314 || mc == 5324 || mc == 5332 || mc == 5334 || mc == 5142 || mc == 5242 || mc == 5412 || mc == 5422 || mc == 5414 || mc == 5424 || mc == 5342 || mc == 5432 || mc == 5434 || mc == 5442 || mc == 5444 || mc == 5512 || mc == 5522 || mc == 5514 || mc == 5524 || mc == 5532 || mc == 5534 || mc == 5542 || mc == 5544 || mc == 5544 ) {return 1;}
    else {return 0;}
    
}

void findCombine(int index, vector<int> &a, vector<int>& tmp, const int& n, vector<vector<int> >& res){
	if(tmp.size() == n/2){
		res.push_back(tmp);
		return;
	}
	for(int i=index; i<=n; i++){
		tmp.push_back(a[i-1]);
		findCombine(i+1, a, tmp, n, res);
		tmp.pop_back();
	}
	return;
}
vector<vector<int> > Combine(int n, vector<int> a){
	vector<vector<int> > res;
	if(n<1 || n%2!=0 || a.size()<1)
		return res;
	vector<int> tmp;
	findCombine(1, a, tmp, n, res);
	return res;
}
vector<int> DifferenceSet(vector<int> total, vector<int> original, vector<int> left){
    for(int i = 0; i<total.size(); i++)
    {
        int element = total[i];
        vector<int >::iterator it;
        it = find(original.begin(), original.end(), element);
        if(it == original.end()){
            left.push_back(element);
        }
    }
    return left;
}
std::vector<vector<int> > pair4jets(int numjets);

// from large to small
static bool sortEn(ReconstructedParticle* a1, ReconstructedParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

// from large to small
static bool sortEnMC(MCParticle* a1, MCParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

void vcb2::init() {
    
    printParameters();
    TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));
    if (!tree_file->IsOpen()) {
        delete tree_file;
        tree_file=new TFile(_treeFileName.c_str(),"NEW");
    }
    
    _outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
    _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
    
    _outputTree->Branch("EventNr",       &eventNr,        "EventNr/I");
    _outputTree->Branch("Num",           &Num,            "Num/I");


    _outputTree->Branch("HbbL",          &HbbL,                "HbbL[2]/F");
    _outputTree->Branch("HccL",          &HccL,                "HccL[2]/F");
    
    _outputTree->Branch("jet14m",        &jet14m,              "jet14m[4]/F");
    _outputTree->Branch("jet24m",        &jet24m,              "jet24m[4]/F");

    _outputTree->Branch("quark14m",      &quark14m,            "quark14m[4]/F");
    _outputTree->Branch("quark24m",      &quark24m,            "quark24m[4]/F");
    
    _outputTree->Branch("quark1PDG",     &quark1PDG,           "quark1PDG/I");
    _outputTree->Branch("quark2PDG",     &quark2PDG,           "quark2PDG/I");
    
//    _outputTree->Branch("jet1LeadHadEnRatio",     &jet1LeadHadEnRatio,           "jet1LeadHadEnRatio/F");
//    _outputTree->Branch("jet2LeadHadEnRatio",     &jet2LeadHadEnRatio,           "jet2LeadHadEnRatio/F");
    
    _outputTree->Branch("lead1HadEn",     &lead1HadEn,           "lead1HadEn/F");
    _outputTree->Branch("lead2HadEn",     &lead2HadEn,           "lead2HadEn/F");
    _outputTree->Branch("lead1PID",     &lead1PID,           "lead1PID/I");
    _outputTree->Branch("lead2PID",     &lead2PID,           "lead2PID/I");
    _outputTree->Branch("lead1EndP",        &lead1EndP,              "lead1EndP[3]/F");
    _outputTree->Branch("lead2EndP",        &lead2EndP,              "lead2EndP[3]/F");
    _outputTree->Branch("num_92BDaus",     &num_92BDaus,           "num_92BDaus/I");

    
    _outputTree->Branch("thrust",        &thrust,              "thrust/D");
    _outputTree->Branch("Y12",           &Y12,                 "Y12/D");
    _outputTree->Branch("Y23",           &Y23,                 "Y23/D");
    _outputTree->Branch("Y34",           &Y34,                 "Y34/D");
 
    Num = 0;
}


bool isQuark(int pdg) {
    return abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6;
}

void getYs(LCEvent * evtP, int &Y12, int &Y23, int &Y34)
{
    Y12 = 999, Y23 = 999, Y34 = 999;

    LCCollection* col_Jet = evtP->getCollection( "RefinedJets" );
    int num_jet = col_Jet->getNumberOfElements();
    cout<<"num_jet : "<<num_jet<<endl;
    
    if(num_jet == 2) {
        PIDHandler pidh(col_Jet);
        int algo   = 0;
        algo = pidh.getAlgorithmID("yth");
        int iy12 = -1, iy23 = -1, iy34 = -1;
        iy12 = pidh.getParameterIndex (algo, "y12" );
        iy23 = pidh.getParameterIndex (algo, "y23" );
        iy34 = pidh.getParameterIndex (algo, "y34" );
        ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
        const ParticleID &pid = pidh.getParticleID(jet, algo);
        Y12 = pid.getParameters()[iy12];
        Y23 = pid.getParameters()[iy23];
        Y34 = pid.getParameters()[iy34];
    }    
}

void vcb2::processEvent( LCEvent * evtP )
{
    
    
    if (evtP)
    {
        try{
            
            cout<<"Next Event *******************************************************************************************************"<<endl;
            eventNr = evtP->getEventNumber();
            cout<<"eventNr : "<<eventNr<<" Num : "<<Num<<endl;
            
            
            //ArborPFOs
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int nPFO = col_PFO->getNumberOfElements();
            cout<<"nPFO : "<<nPFO<<endl;
            std::vector<TLorentzVector>         vArborTL;   vArborTL.clear();

            for(int i = 0; i<nPFO; i++){
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                vArborTL.push_back(temp);
            }

            // quark list ...
            std::vector<MCParticle*> quarkvec;
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int const n_MCP = col_MCP->getNumberOfElements();
            
            std::vector<MCParticle*> v92Daus;

//            std::map<MCParticle*, float> mapMCcount; mapMCcount.clear();
            
            for(int i = 0; i < n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int const NParents      = a_MCP->getParents().size();
                int const NDaughters    = a_MCP->getDaughters().size();
                int const PDG           = a_MCP->getPDG();
                

                
                if(NParents == 0 && isQuark(PDG) ) {
                    quarkvec.push_back(a_MCP);
                }
                
                if(PDG == 92) {
                    for(int j = 0; j < NDaughters; j++){
                        MCParticle* mcp_92Dau = a_MCP->getDaughters()[j];
                        int tempPDG = mcp_92Dau->getPDG();
                        if(ISB(tempPDG) == 1 || ISBbar(tempPDG) == 1) {                            
                            if(v92Daus.size() == 0) { v92Daus.push_back(mcp_92Dau); }
                            else{
                                int count = 0;
                                for(int k = 0; k<v92Daus.size(); k++){
                                    MCParticle* temp = v92Daus.at(k);
                                    if( abs(temp->getEnergy() - mcp_92Dau->getEnergy())<0.01 ){ count = 1; }
                                }
                                if( count == 0 ){ v92Daus.push_back(mcp_92Dau); }
                            }
                            
                            
                        }
                    }
                }
                
                
                
            }
                        
            for(int i = 0; i<4; i++){
                jet14m[i] = 999;
                jet24m[i] = 999;
                
                quark14m[i] = 999;
                quark24m[i] = 999;
            }
            
            
            sort(v92Daus.begin(), v92Daus.end(), sortEnMC);
            
            if(v92Daus.size() != 0){
                for(int i = 0; i<v92Daus.size(); i++){
                    MCParticle* a_MCP = v92Daus.at(i);
                    cout<<"a_MCP->getPDG() : "<<a_MCP->getPDG()<<endl;
                }
            
            }
            
            
            MCParticle* lead1 = NULL; MCParticle* lead2 = NULL;
            MCParticle* lead1LeadHad = NULL; MCParticle* lead2LeadHad = NULL;
            num_92BDaus = 0; num_92BDaus = v92Daus.size();
            cout<<"num_92BDaus : "<<num_92BDaus<<endl;
            
            if(v92Daus.size() > 1){
                lead1 = v92Daus.at(0);
                lead2 = v92Daus.at(1);
                int stop = 1, immediately = 0;
                
                while( lead1->getDaughters().size() != 0 && stop != 0){
                    stop = 0;
//                    immediately += 1;
                    for(int i = 0; i<lead1->getDaughters().size(); i++){
                        MCParticle* Dau = lead1->getDaughters()[i];
                        int DauPDG = Dau->getPDG();
                        if(ISB(DauPDG) == 1 || ISBbar(DauPDG) == 1){ lead1 = Dau; stop += 1;}
                    }
                }
//                if(stop == 0 ){ lead1LeadHad = lead1; }
                
                
                stop = 1; immediately = 0;
                while( lead2->getDaughters().size() != 0 && stop != 0){
                    stop = 0;
//                    immediately += 1;
                    for(int i = 0; i<lead2->getDaughters().size(); i++){
                        MCParticle* Dau = lead2->getDaughters()[i];
                        int DauPDG = Dau->getPDG();
                        if(ISB(DauPDG) == 1 || ISBbar(DauPDG) == 1){ lead2 = Dau; stop += 1;}
                    }
                }
                

                cout<<"lead1LeadHad->getPDG() : "<<lead1->getPDG()<<" : "<<lead2->getPDG()<<endl;
                            
            }
                
            
            
            quark1PDG = 999, quark2PDG = 999;
            if(quarkvec.size() == 2){
                MCParticle* quark1 = quarkvec.at(0);
                MCParticle* quark2 = quarkvec.at(1);
                quark1PDG = quark1->getPDG();
                quark2PDG = quark2->getPDG();

                quark14m[0] = quark1->getMomentum()[0];
                quark14m[1] = quark1->getMomentum()[1];
                quark14m[2] = quark1->getMomentum()[2];
                quark14m[3] = quark1->getEnergy();
                
                quark24m[0] = quark2->getMomentum()[0];
                quark24m[1] = quark2->getMomentum()[1];
                quark24m[2] = quark2->getMomentum()[2];
                quark24m[3] = quark2->getEnergy();
                
            }
            
            
            getYs(evtP, Y12, Y23, Y34);
            
            for(int i = 0; i<2; i++){
                HbbL[i] = 999;
                HccL[i] = 999;
            }
            

            
//            jet1LeadHadEnRatio = 0, jet2LeadHadEnRatio = 0;
            lead1HadEn = 0, lead2HadEn = 0;
            
            for(int i = 0; i < 3; i++){
                lead1EndP[i] = 999;
                lead2EndP[i] = 999;
            }

            lead1PID = 999, lead2PID = 999;
  
            if(num_jet == 2){
                
                PIDHandler pidh(col_Jet);
                int algo   = 0;
                int ibtag  =-1;
                int ictag  =-1;
                int ibctag =-1;
                int icat   =-1;
                algo      = pidh.getAlgorithmID("lcfiplus");
                ibtag     = pidh.getParameterIndex (algo,  "BTag");
                ictag     = pidh.getParameterIndex (algo,  "CTag");
                ibctag    = pidh.getParameterIndex (algo,  "BCTag");
                icat      = pidh.getParameterIndex (algo,  "Category");
                
                ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
                TLorentzVector TLjet1(jet1->getMomentum(), jet1->getEnergy());
                TVector3 TVjet1 = TLjet1.Vect();
                
                const ParticleID &pid1 = pidh.getParticleID(jet1, algo);
                double btag1  = pid1.getParameters()[ibtag];
                double ctag1  = pid1.getParameters()[ictag];

                ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(1));
                TLorentzVector TLjet2(jet2->getMomentum(), jet2->getEnergy());
                TVector3 TVjet2 = TLjet2.Vect();
                
                
                if( v92Daus.size() > 1 ){
//                    lead1 = lead1LeadHad;
//                    lead2 = lead2LeadHad;
                    
                    TLorentzVector TLlead1(lead1->getMomentum(), lead1->getEnergy());
                    TVector3 TVlead1 = TLlead1.Vect();
                    TLorentzVector TLlead2(lead2->getMomentum(), lead2->getEnergy());
                    TVector3 TVlead2 = TLlead2.Vect();
                    
                    if(TVlead1.Angle(TVjet1) + TVlead2.Angle(TVjet2) < TVlead1.Angle(TVjet2) + TVlead2.Angle(TVjet1) ){
                        lead1HadEn = lead1->getEnergy(); lead1PID = lead1->getPDG();
                        TVector3 lead1EP = lead1->getEndpoint();
                        lead1EndP[0] = lead1EP[0]; lead1EndP[1] = lead1EP[1]; lead1EndP[2] = lead1EP[2];
                        
                        lead2HadEn = lead2->getEnergy(); lead2PID = lead2->getPDG();
                        TVector3 lead2EP = lead2->getEndpoint();
                        lead2EndP[0] = lead2EP[0]; lead2EndP[1] = lead2EP[1]; lead2EndP[2] = lead2EP[2];

                    }
                    else{
                        TVector3 lead1EP = lead1->getEndpoint();
                        TVector3 lead2EP = lead2->getEndpoint();

                        lead1HadEn = lead2->getEnergy(); lead1PID = lead2->getPDG();
                        lead1EndP[0] = lead2EP[0]; lead1EndP[1] = lead2EP[1]; lead1EndP[2] = lead2EP[2];

                        lead2HadEn = lead1->getEnergy(); lead2PID = lead1->getPDG();
                        lead2EndP[0] = lead1EP[0]; lead2EndP[1] = lead1EP[1]; lead2EndP[2] = lead1EP[2];

                    }
                
                }
                
                
                const ParticleID &pid2 = pidh.getParticleID(jet2, algo);
                double btag2  = pid2.getParameters()[ibtag];
                double ctag2  = pid2.getParameters()[ictag];

                HbbL[0] = btag1; HbbL[1] = btag2;
                HccL[0] = ctag1; HccL[1] = ctag2;
                
                jet14m[0] = jet1->getMomentum()[0];
                jet14m[1] = jet1->getMomentum()[1];
                jet14m[2] = jet1->getMomentum()[2];
                jet14m[3] = jet1->getEnergy();
                
                jet24m[0] = jet2->getMomentum()[0];
                jet24m[1] = jet2->getMomentum()[1];
                jet24m[2] = jet2->getMomentum()[2];
                jet24m[3] = jet2->getEnergy();
                
                
                
                if(quarkvec.size() == 2){
                    MCParticle* quark1 = quarkvec.at(0);
                    MCParticle* quark2 = quarkvec.at(1);
                    
                    TLorentzVector TLquark1(quark1->getMomentum(), quark1->getEnergy());
                    TLorentzVector TLquark2(quark2->getMomentum(), quark2->getEnergy());
                    TVector3 TVquark1 = TLquark1.Vect();
                    TVector3 TVquark2 = TLquark2.Vect();
                    
                    if(TVjet1.Angle(TVquark1) + TVjet2.Angle(TVquark2) < TVjet1.Angle(TVquark2) + TVjet2.Angle(TVquark1)){
                        
                        quark1PDG = quark1->getPDG();
                        quark2PDG = quark2->getPDG();
                        
                        quark14m[0] = quark1->getMomentum()[0]; quark14m[1] = quark1->getMomentum()[1];
                        quark14m[2] = quark1->getMomentum()[2]; quark14m[3] = quark1->getEnergy();
                        quark24m[0] = quark2->getMomentum()[0]; quark24m[1] = quark2->getMomentum()[1];
                        quark24m[2] = quark2->getMomentum()[2]; quark24m[3] = quark2->getEnergy();
                    }
                    else{
                        
                        quark1PDG = quark2->getPDG();
                        quark2PDG = quark1->getPDG();
                        
                        quark14m[0] = quark2->getMomentum()[0]; quark14m[1] = quark2->getMomentum()[1];
                        quark14m[2] = quark2->getMomentum()[2]; quark14m[3] = quark2->getEnergy();
                        quark24m[0] = quark1->getMomentum()[0]; quark24m[1] = quark1->getMomentum()[1];
                        quark24m[2] = quark1->getMomentum()[2]; quark24m[3] = quark1->getEnergy();
                    }
                    
                }
                
            }
            
            
            for(int i = 0; i<4; i++){
                cout<<"jet14m[i] : "<<jet14m[i]<<endl;
                cout<<"jet24m[i] : "<<jet24m[i]<<endl;

                cout<<"quark14m[i] : "<<quark14m[i]<<endl;
                cout<<"quark24m[i] : "<<quark24m[i]<<endl;

            }
            
            for(int i = 0; i<3; i++){
                cout<<"lead1EndP[i] : "<<lead1EndP[i]<<endl;
                cout<<"lead2EndP[i] : "<<lead2EndP[i]<<endl;
            }
            
            cout<<"lead1PID : "<< lead1PID <<" lead2PID : "<<lead2PID<<endl;
            

        }catch (lcio::DataNotAvailableException err) {  }

    }
    
    _outputTree->Fill();
    Num ++;
}



void vcb2::end()
{
    
    if (_outputTree) {
        
        TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        //tree_file->cd();
        tree_file->Write();
        delete tree_file;
        //tree_file->Close();
    }
    
}


