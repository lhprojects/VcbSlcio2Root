#include <vcb.hh>
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

vcb a_vcb_instance;


vcb::vcb()
: Processor("vcb"),
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
    
    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             "MCPSIMUFSP",
                             "the final state particle after simulation",
                             _outmcpsimufsp,
                             std::string("ReconstructedParticle") );
    
    
    registerOutputCollection( LCIO::LCRELATION,
                             "mcpsimurelation",
                             " relation between MCP and Reco after simulation",
                             _outMCPSIMURelation,
                             std::string("RelationusedSed") );
    
    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}





static bool sortEn(ReconstructedParticle* a1, ReconstructedParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

void vcb::init() {
    
    printParameters();
//    TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));
//    if (!tree_file->IsOpen()) {
//        delete tree_file;
//        tree_file=new TFile(_treeFileName.c_str(),"NEW");
//    }
    
//    _outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
//    _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
//    _outputTree->Branch("EventNr",     &eventNr,       "EventNr/I");
//    _outputTree->Branch("Num",         &Num,           "Num/I");
//    _outputTree->Branch("visEn",        &visEn,         "visEn/D");
//    _outputTree->Branch("multiplicity", &multiplicity,  "multiplicity/I");
//    _outputTree->Branch("LeadElecEn",   &LeadElecEn,    "LeadElecEn/D");
//    _outputTree->Branch("leadMuonEn",   &leadMuonEn,    "leadMuonEn/D");
//    _outputTree->Branch("subleadMuonEn",   &subleadMuonEn,    "subleadMuonEn/D");
//    _outputTree->Branch("ratio",    &ratio,    "ratio/D");


 
    Num = 0;
}

void vcb::processEvent( LCEvent * evtP )
{
    
    
    LCCollectionVec* otmcpsimufsp = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
    LCCollectionVec* otMCPSIMURelation = new LCCollectionVec( LCIO::LCRELATION );
    LCFlagImpl newMCPSIMUlinkflag;
    newMCPSIMUlinkflag.setBit(LCIO::CHBIT_LONG);
    otMCPSIMURelation->setFlag(newMCPSIMUlinkflag.getFlag());
    
    if (evtP)
    {
        try{
            
            cout<<"Next Event *******************************************************************************************************"<<endl;
//            eventNr = evtP->getEventNumber();
//            cout<<"eventNr : "<<eventNr<<" Num : "<<Num<<endl;
            
            //ArborPFOs
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int nPFO = col_PFO->getNumberOfElements();
            cout<<"nPFO : "<<nPFO<<endl;
            TLorentzVector TLPFO(0,0,0,0);
            
            
            std::vector<ReconstructedParticle*> vElec;      vElec.clear();
            std::vector<ReconstructedParticle*> vMuon;      vMuon.clear();

            for(int i = 0; i<nPFO; i++){
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                TVector3 TVtemp = temp.Vect();
                TLPFO += temp;
                int type = abs(pfo->getType());
                if(type == 11){vElec.push_back(pfo);}
                if(type == 13){vMuon.push_back(pfo);}
                else{
                    ReconstructedParticleImpl* a_Reco = new ReconstructedParticleImpl();
                    a_Reco->addParticle(pfo);
                    a_Reco->setEnergy( pfo->getEnergy() );
                    a_Reco->setMomentum( pfo->getMomentum() );
                    a_Reco->setType( pfo->getType() );
                    if( pfo->getTracks().size() ){
                        a_Reco->addTrack(pfo->getTracks()[0]);
                    }
                    
                    a_Reco->setCharge(pfo->getCharge());
                    a_Reco->setCovMatrix(pfo->getCovMatrix());
                    if(pfo->getClusters().size()){
                        for(int j = 0; j<pfo->getClusters().size(); j++){
                            a_Reco->addCluster( pfo->getClusters()[j] );
                        }
                    }
                    otmcpsimufsp->addElement( a_Reco );
                    LCRelationImpl *newrel = new LCRelationImpl(a_Reco, pfo, 1.0);
                    otMCPSIMURelation->addElement( newrel );
                }
                
                
            }
            
            cout<<"the total energy of ArborPFOs is : "<<TLPFO.E()<<endl;
//            TVector3 TVPFO = TLPFO.Vect();
//            Pt = TVPFO.Perp();
//            Pl = pow(TVPFO.Mag2() - TVPFO.Perp2(), 0.5);


            
            
           
                
            //the following code used to analysis the multiplicity information of final state particles
//            visEn = 0; multiplicity = 0;
//            visEn = TLPFO.E();
//            multiplicity = nPFO;
            
//            double leadMuonEn = 0, LeadElecEn = 0, subleadMuonEn = 0, ratio = 0;
            

            
//            sort(vElec.begin(), vElec.end(), sortEn);
//            sort(vMuon.begin(), vMuon.end(), sortEn);
//            
//
//            ReconstructedParticle* leadMuon = NULL;
//            ReconstructedParticle* leadElec = NULL;
//            if(vMuon.size() > 0){
//                leadMuon = vMuon.at(0);
// //               leadMuonEn = leadMuon->getEnergy();
//                
//                TLorentzVector TLseedMuon(leadMuon->getMomentum(), leadMuon->getEnergy());
//                TVector3 TVseedMuon = TLseedMuon.Vect();
//                
//                TLorentzVector TLaround(0, 0, 0, 0);
//                for(int i = 0; i<nPFO; i++){
//                    ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
//                    TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
//                    TVector3 TVtemp = temp.Vect();
//                    
//                    if(pfo != leadMuon){
//                        if( TVtemp.Angle(TVseedMuon) < 0.5 ){ TLaround += temp; }
//                    }
//                    
//                }
//                
//                ratio = leadMuon->getEnergy()/( TLaround.E() + leadMuon->getEnergy() );
//                
//                if(vMuon.size() > 1){
//                    ReconstructedParticle* subleadMuon = vMuon.at(1);
//                    subleadMuonEn = subleadMuon->getEnergy();
//                }
//            }
//            
//            if(vElec.size() > 0){
//                leadElec = vElec.at(0);
//                LeadElecEn = leadElec->getEnergy();
//            }
//            
//       
            
            if(vMuon.size() > 1){
                for(int i = 1; i<vMuon.size(); i++){
                    ReconstructedParticle* pfo = vMuon.at(i);
                    ReconstructedParticleImpl* a_Reco = new ReconstructedParticleImpl();
                    a_Reco->addParticle(pfo);
                    a_Reco->setEnergy( pfo->getEnergy() );
                    a_Reco->setMomentum( pfo->getMomentum() );
                    a_Reco->setType( pfo->getType() );
                    if( pfo->getTracks().size() ){
                        a_Reco->addTrack(pfo->getTracks()[0]);
                    }
                    
                    a_Reco->setCharge(pfo->getCharge());
                    a_Reco->setCovMatrix(pfo->getCovMatrix());
                    if(pfo->getClusters().size()){
                        for(int j = 0; j<pfo->getClusters().size(); j++){
                            a_Reco->addCluster( pfo->getClusters()[j] );
                        }
                    }
                    otmcpsimufsp->addElement( a_Reco );
                    LCRelationImpl *newrel = new LCRelationImpl(a_Reco, pfo, 1.0);
                    otMCPSIMURelation->addElement( newrel );
                }
            }
            
            cout<<"nPFO : "<<nPFO<<" vMuon.size() : "<<vMuon.size()<<endl;
            cout<<"otmcpsimufsp->getNumberOfElements() : "<<otmcpsimufsp->getNumberOfElements()<<endl;
            cout<<"otMCPSIMURelation->getNumberOfElements() : "<<otMCPSIMURelation->getNumberOfElements()<<endl;
            

        }catch (lcio::DataNotAvailableException err) {  }

    }
    
    evtP->addCollection( otmcpsimufsp, _outmcpsimufsp.c_str() );
    evtP->addCollection( otMCPSIMURelation, _outMCPSIMURelation.c_str() );
    
//    _outputTree->Fill();
    Num ++;
}



void vcb::end()
{
    
    if (_outputTree) {
        
        TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        //tree_file->cd();
        tree_file->Write();
        delete tree_file;
        //tree_file->Close();
    }
    
}





