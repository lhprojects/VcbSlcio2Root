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

int ISBbar(int mc){
    if(mc == 511 || mc == 521 || mc == 10511 || mc == 10521 || mc == 513 || mc == 523 || mc == 10513 || mc == 10523 || mc == 20513 || mc == 20523 || mc == 515 || mc == 525 || mc == 531 || mc == 10531 || mc == 533 || mc == 10533 || mc == 20533 || mc == 535 || mc == 541 || mc == 10541 || mc == 543 || mc == 10543 || mc == 20543 || mc == 545 || mc == -5122 || mc == -5112 || mc == -5212 || mc == -5222 || mc == -5114 || mc == -5214 || mc == -5224 || mc == -5132 || mc == -5232 || mc == -5312 || mc == -5322 || mc == -5314 || mc == -5324 || mc == -5332 || mc == -5334 || mc == -5142 || mc == -5242 || mc == -5412 || mc == -5422 || mc == -5414 || mc == -5424 || mc == -5342 || mc == -5432 || mc == -5434 || mc == -5442 || mc == -5444 || mc == -5512 || mc == -5522 || mc == -5514 || mc == -5524 || mc == -5532 || mc == -5534 || mc == -5542 || mc == -5544 || mc == -5544 ) {return 1;}
    else {return 0;}
}

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

static bool sortEn(ReconstructedParticle* a1, ReconstructedParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

static bool sortEnMC(MCParticle* a1, MCParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

void CalcuThrust(std::vector<TLorentzVector > UsedForThrust, std::vector<double> &result);


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
            
            thrust = 999;
            std::vector<double> shape1Result; shape1Result.clear();
            CalcuThrust(vArborTL, shape1Result);
            thrust = shape1Result.at(0);

            std::vector<MCParticle*> quarkvec; quarkvec.clear();
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            
            std::vector<MCParticle*> v92Daus;      v92Daus.clear();

//            std::map<MCParticle*, float> mapMCcount; mapMCcount.clear();
            
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents      = a_MCP->getParents().size();
                int NDaughters    = a_MCP->getDaughters().size();
                int PDG           = a_MCP->getPDG();
                

                
                if(NParents == 0 && (abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6) ){
                    quarkvec.push_back(a_MCP);
                }
                
                if(PDG == 92){
                    cout<<"92...."<<endl;
                    for(int j = 0; j<NDaughters; j++){
                        MCParticle* mcp_92Dau = a_MCP->getDaughters()[j];
                        int tempPDG = mcp_92Dau->getPDG();
                        if(ISB(tempPDG) == 1 || ISBbar(tempPDG) == 1){
//                            mapMCcount[mcp_92Dau] = float;
                            cout<<"92 B daughters ..."<<tempPDG<<" mcp_92Dau->getEnergy() : "<<mcp_92Dau->getEnergy()<<endl;
//                            if( std::find(v92Daus.begin(), v92Daus.end(), mcp_92Dau) == v92Daus.end() ){
//                                v92Daus.push_back(mcp_92Dau);
//                            }
                            
                            if(v92Daus.size() == 0){ v92Daus.push_back(mcp_92Dau); }
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
            
            
 //           cout<<"mapMCcount.size() : "<<mapMCcount.size()<<endl;
            
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
//            std::vector<MCParticle*> lead1Daus; lead1Daus.clear();
//            std::vector<MCParticle*> lead2Daus; lead2Daus.clear();
            
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
//                if(stop == 0 ){ lead2LeadHad = lead2; }
                

                cout<<"lead1LeadHad->getPDG() : "<<lead1->getPDG()<<" : "<<lead2->getPDG()<<endl;
                
                
//                for(int i = 0; i<n_MCP; i++)
//                {
//                    MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
//                    int NParents      = a_MCP->getParents().size();
//                    int PDG           = a_MCP->getPDG();
//                    
//                    MCParticle* a_parent = a_MCP;
//                    if(abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16 && abs(PDG) != 11 && abs(PDG) != 13 && abs(PDG) != 15 && abs(PDG) != 22 && a_MCP->getGeneratorStatus() == 1 ){
//                        while( a_parent->getParents().size() != 0 && a_parent != lead1 && a_parent != lead2 ){
//                            a_parent = a_parent->getParents()[0];
//                        }
//                        if( a_parent == lead1 ){lead1Daus.push_back(a_MCP);}
//                        else if(a_parent == lead2){lead2Daus.push_back(a_MCP);}
//                    }
//                }
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
            
            
            
            for(int i = 0; i<2; i++){
                HbbL[i] = 999;
                HccL[i] = 999;
            }
            


            LCCollection* col_Jet = evtP->getCollection( "RefinedJets" );
            int num_jet = col_Jet->getNumberOfElements();
            cout<<"num_jet : "<<num_jet<<endl;
            
            Y12 = 999, Y23 = 999, Y34 = 999;
            if(num_jet == 2){
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
            
//            jet1LeadHadEnRatio = 0, jet2LeadHadEnRatio = 0;
            lead1HadEn = 0, lead2HadEn = 0;
            
            for(int i = 0; i<3; i++){
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
            
            cout<<"lead1PID : "<<lead1PID<<" lead2PID : "<<lead2PID<<endl;
            

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


void CalcuThrust(std::vector<TLorentzVector > UsedForThrust, std::vector<double> &result){
    result.clear();
    double T = 0;
    //the following code used to find the thrust
    double thetaMin = TMath::Pi(), phiMin = 2*TMath::Pi();
    double thetaMin2 = 0, phiMin2 = 0;
    double thetaL = 0, phiL = 0, thetaR = TMath::Pi(), phiR = 2*TMath::Pi();
    int iter = 0;
    double Told = 0;
    double Tnew = 0;
    double cut = 1;
    double thetaRange = 0, phiRange = 0;
    do{
        iter += 1;
        //        cout<<"iter : "<<iter<<endl;
        if(iter == 1){
            thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
        }
        else if(iter != 1){
            
            thetaRange = 0.1*(thetaR - thetaL);
            phiRange = 0.1*(phiR - phiL);
            
            thetaL =  thetaMin - thetaRange;
            thetaR = thetaMin + thetaRange;
            phiL = phiMin - phiRange;
            phiR = phiMin + phiRange;
            thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
            
            //            cout<<"thetaL : "<<thetaL<<" thetaR : "<<thetaR<<endl;
            //            cout<<"phiL : "<<phiL<<" phiR : "<<phiR<<endl;
        }
        //        cout<<"thetaRange : "<<thetaRange<<" phiRange : "<<phiRange<<endl;
        for(double theta = thetaL; theta <= thetaR; theta += 0.1*thetaRange){   //in this round, find the max T
            for(double phi = phiL; phi <= phiR; phi += 0.1*phiRange){
                
                double x = sin(theta)*cos(phi);
                double y = sin(theta)*sin(phi);
                double z = cos(theta);
                
                double denominator = 0;
                double numerator = 0;
                for(int i = 0; i<UsedForThrust.size(); i++){
                    TLorentzVector TLtemp = UsedForThrust.at(i);
                    //                    TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                    TVector3 TVtemp = TLtemp.Vect();
                    denominator += TVtemp.Mag();
                    numerator += abs(x*TVtemp(0) + y*TVtemp(1) + z*TVtemp(2));
                }
                double Ttemp = numerator/denominator;
                if(Ttemp > T){
                    thetaMin = theta;   phiMin = phi; T = Ttemp;
                    //                   cout<<"*************"<<endl;
                    //                   cout<<"T : "<<T<<"thetaMin : phiMin "<<thetaMin<<" : "<<phiMin<<endl;
                    //                   cout<<"*************"<<endl;
                }
            }
        }
        if(iter == 1){Told = T; Tnew = T;}
        else if(T >= Tnew && iter != 1){
            Told = Tnew; Tnew = T; cut = (Tnew - Told)/Tnew;
        }
        //        cout<<"cut : "<<cut<<endl;
    }
    while(cut >= 0.2);
    
    //    result[3] = {T, phiMin, thetaMin};
    result.push_back(T);
    
    TVector3 tempThrust(0,0,0);
    tempThrust.SetXYZ(sin(thetaMin)*cos(phiMin), sin(thetaMin)*sin(phiMin), cos(thetaMin));
    
    //the following code used to get Hemisphere masses
    std::vector<TLorentzVector > hemisphere1;
    std::vector<TLorentzVector > hemisphere2;
    hemisphere1.clear(); hemisphere2.clear();
    double visEn = 0;
    double JetBroadeningDenominator = 0;
    TVector3 TVtotal(0,0,0);
    for(int i = 0; i<UsedForThrust.size(); i++){
        TLorentzVector TLtemp = UsedForThrust.at(i);
        //        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TVtotal += TVtemp;
        if(TVtemp.Angle(tempThrust) > 0.5*TMath::Pi()){hemisphere1.push_back(TLtemp);}
        else {hemisphere2.push_back(TLtemp);}
        visEn += TLtemp.E();
        JetBroadeningDenominator += TVtemp.Mag();
    }
    double hemi1En = 0, hemi1Mass = 0, hemi2En = 0, hemi2Mass = 0;
    TLorentzVector TLsphere1(0,0,0,0);
    TLorentzVector TLsphere2(0,0,0,0);
    for(int i = 0; i<hemisphere1.size(); i++){
        TLsphere1 += hemisphere1.at(i);
    }
    for(int i = 0; i<hemisphere2.size(); i++){
        TLsphere2 += hemisphere2.at(i);
    }
    hemi1En = TLsphere1.E();   hemi1Mass = TLsphere1.M();
    hemi2En = TLsphere2.E();   hemi2Mass = TLsphere2.M();
    cout<<"hemi1En : "<<hemi1En<<" hemi2En : "<<hemi2En<<endl;
    cout<<"hemi1Mass : "<<hemi1Mass<<" hemi2Mass : "<<hemi2Mass<<endl;
    //    cout<<"the number of particles in two hemispheres is "<<hemisphere1.size()+hemisphere2.size()<<endl;
    
    
    double JetBroadeningNumerator1 = 0, JetBroadeningNumerator2 =0;
    TLorentzVector TLHemi1(0,0,0,0);
    TLorentzVector TLHemi2(0,0,0,0);
    for(int i = 0; i<hemisphere1.size(); i++){
        TLorentzVector TLtemp = hemisphere1.at(i);
        //        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TLHemi1 += TLtemp;
        JetBroadeningNumerator1 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
        
    }
    for(int i = 0; i<hemisphere2.size(); i++){
        TLorentzVector TLtemp = hemisphere2.at(i);
        //        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TLHemi2 += TLtemp;
        JetBroadeningNumerator2 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
    }
    double hemiMass1 = 0, hemiMass2 = 0, hemiBroadening1 = 0, hemiBroadening2 = 0;
    if(visEn != 0){
        hemiMass1 = (TLHemi1.M())*(TLHemi1.M())/(visEn*visEn);
        hemiMass2 = (TLHemi2.M())*(TLHemi2.M())/(visEn*visEn);
    }
    if(JetBroadeningDenominator != 0){
        hemiBroadening1 = JetBroadeningNumerator1/(2*JetBroadeningDenominator);
        hemiBroadening2 = JetBroadeningNumerator2/(2*JetBroadeningDenominator);
    }
    
    //    cout<<"hemiMass1 : "<<hemiMass1<<" hemiMass2 : "<<hemiMass2<<endl;
    //    cout<<"hemiBroadening1 : "<<hemiBroadening1<<" hemiBroadening2 : "<<hemiBroadening2<<endl;
    result.push_back(hemiMass1);
    result.push_back(hemiMass2);
    result.push_back(hemiBroadening1);
    result.push_back(hemiBroadening2);
    
    
    
    //the following code used to get rapidity
    TVector3 TVParaThrust = tempThrust.Orthogonal();
    double transverseM = TVtotal.Perp(TVParaThrust);
    double rapidity = 0.5*log((visEn + transverseM)/(visEn - transverseM));
    //    result.push_back(rapidity);
    result.push_back( tempThrust.CosTheta() );
    result.push_back( hemi1En ); result.push_back(hemi2En); result.push_back(hemi1Mass); result.push_back(hemi2Mass);
    
    cout<<"visEn : "<<visEn<<endl;
    
    //    return 0;
}




