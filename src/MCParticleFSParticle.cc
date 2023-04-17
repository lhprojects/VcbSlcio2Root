
#include "MCParticleFSParticle.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>
#include <algorithm>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <marlin/VerbosityLevels.h>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

std::vector<MCParticle*> getResult(std::vector<MCParticle*> vec, std::vector<MCParticle*> vec92);

MCParticleFSParticle aMCParticleFSParticle ;

MCParticleFSParticle::MCParticleFSParticle()
	: Processor("MCParticleFSParticle") {

		_description = "MCP final state particle finding" ;


        registerInputCollection( LCIO::MCPARTICLE,
                                "InputCollection" ,
                                "Input collection of MCparticles",
                                _inputPFOsCollection,
                                std::string("MCParticle"));
        
        registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                                 "MCPSIMUFSP",
                                 "the final state particle after simulation",
                                 _outmcpsimufsp,
                                 std::string("GenJetParticles") );
        

        registerOutputCollection( LCIO::LCRELATION,
                                 "mcpsimurelation",
                                 " relation between MCP and Reco after simulation",
                                 _outMCPSIMURelation,
                                 std::string("RelationusedSed") );
        
        registerProcessorParameter( "replaceFL",
                        "replace final state particles come from PDG of replaceFL",
                        _replaceFL,
                        int(9999999999));
	}


void MCParticleFSParticle::init() {
	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	printParameters() ;
}

void MCParticleFSParticle::processRunHeader( LCRunHeader* run) {
}

static bool sortEn(MCParticle* a1, MCParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

bool cmp(const pair<MCParticle*, double>& a, const pair<MCParticle*, double>& b){
    return a.second < b.second;
}

void MCParticleFSParticle::processEvent( LCEvent * evt ) {


    std::cout << "*************************************" << std::endl;  
    std::cout << "MCPFinalStateParticle: processEvent " << std::endl;
    std::vector<int> allkinds92Dau; allkinds92Dau.clear();
    allkinds92Dau.push_back(99999); allkinds92Dau.push_back(213);   allkinds92Dau.push_back(10213); allkinds92Dau.push_back(113);
    allkinds92Dau.push_back(223);   allkinds92Dau.push_back(211);   allkinds92Dau.push_back(10323); allkinds92Dau.push_back(10313); allkinds92Dau.push_back(321);
    allkinds92Dau.push_back(215);   allkinds92Dau.push_back(311);   allkinds92Dau.push_back(323);   allkinds92Dau.push_back(313);   allkinds92Dau.push_back(10113);
    allkinds92Dau.push_back(10223); allkinds92Dau.push_back(2212);  allkinds92Dau.push_back(2112);  allkinds92Dau.push_back(111);   allkinds92Dau.push_back(221);
    allkinds92Dau.push_back(115);   allkinds92Dau.push_back(225);   allkinds92Dau.push_back(20213); allkinds92Dau.push_back(10211); allkinds92Dau.push_back(3122);
    allkinds92Dau.push_back(325);   allkinds92Dau.push_back(315);   allkinds92Dau.push_back(331);   allkinds92Dau.push_back(20113); allkinds92Dau.push_back(20223);
    allkinds92Dau.push_back(10333); allkinds92Dau.push_back(10221); allkinds92Dau.push_back(10111); allkinds92Dau.push_back(333);   allkinds92Dau.push_back(2214);
    allkinds92Dau.push_back(1114);  allkinds92Dau.push_back(2224);  allkinds92Dau.push_back(2114);  allkinds92Dau.push_back(20323); allkinds92Dau.push_back(3222);
    allkinds92Dau.push_back(10311); allkinds92Dau.push_back(10321); allkinds92Dau.push_back(3212);  allkinds92Dau.push_back(3112);  allkinds92Dau.push_back(20313);
    allkinds92Dau.push_back(3312);  allkinds92Dau.push_back(3322);  allkinds92Dau.push_back(335);   allkinds92Dau.push_back(3114);  allkinds92Dau.push_back(3214);
    allkinds92Dau.push_back(3224);  allkinds92Dau.push_back(10331); allkinds92Dau.push_back(20333); allkinds92Dau.push_back(3324);  allkinds92Dau.push_back(3314);
    allkinds92Dau.push_back(3334);  allkinds92Dau.push_back(413);   allkinds92Dau.push_back(411);   allkinds92Dau.push_back(423);   allkinds92Dau.push_back(20533);
    allkinds92Dau.push_back(513);

	_pfoCol = evt->getCollection( _inputPFOsCollection ) ;
    
    LCCollectionVec* otmcpsimufsp = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
    LCCollectionVec* otMCPSIMURelation = new LCCollectionVec( LCIO::LCRELATION );
    LCFlagImpl newMCPSIMUlinkflag;
    newMCPSIMUlinkflag.setBit(LCIO::CHBIT_LONG);
    otMCPSIMURelation->setFlag(newMCPSIMUlinkflag.getFlag());

//    std::vector<int> exclude; exclude.clear();
    
    int npfo = _pfoCol->getNumberOfElements();
    int num92 = 0, countISR = 0;
    int numFL = 0, count94 = 0;

    std::vector<MCParticle*> vISR;       vISR.clear();
    std::vector<MCParticle*> v92;
    TLorentzVector TLISR(0,0,0,0);
    TLorentzVector TL92(0,0,0,0);
    TLorentzVector TLVisble(0,0,0,0);
    TLorentzVector TLNeutrino(0,0,0,0);
    for(int i = 0; i<npfo; i++){
        MCParticle* a_MCP = dynamic_cast<MCParticle* >( _pfoCol->getElementAt(i));
        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        int PDG = abs(a_MCP->getPDG());
        if(PDG == 22 && a_MCP->getParents().size() == 0){  TLISR += TLtemp; vISR.push_back(a_MCP);}
        if(PDG == 92){
            TL92 += TLtemp;
            num92 += 1;
            v92.push_back(a_MCP);
        }
        if( a_MCP->getDaughters().size() == 0 ){
            if( abs(a_MCP->getPDG()) != 12 && abs(a_MCP->getPDG()) != 14 && abs(a_MCP->getPDG()) != 16 ){
                TLVisble += TLtemp;
                numFL += 1;
            }
            else{ TLNeutrino += TLtemp; }
        }
        
        if(a_MCP->getDaughters().size() == 0 && a_MCP->getParents().size() != 0){
            MCParticle* a_parent = a_MCP;
            do{ a_parent = a_parent->getParents()[0]; }
            while( a_parent->getParents().size() != 0 && a_parent->getPDG() != 92 && abs(a_parent->getPDG()) != 1 && abs(a_parent->getPDG()) != 2 && abs(a_parent->getPDG()) != 3 && abs(a_parent->getPDG()) !=4 && abs(a_parent->getPDG()) != 5 && abs(a_parent->getPDG()) != 6 && abs(a_parent->getPDG()) != 94 );
            if( abs(a_parent->getPDG()) == 94 || abs(a_parent->getPDG()) == 1 || abs(a_parent->getPDG()) == 2 || abs(a_parent->getPDG()) == 3 || abs(a_parent->getPDG()) ==4 || abs(a_parent->getPDG()) == 5 || abs(a_parent->getPDG()) == 6){count94 = 1;}
        }

        
    }
  
    cout<<"(TLISR + TL92).E() : "<<(TLISR + TL92).E()<<endl;
    cout<<"(TLNeutrino + TLVisble).E() : "<<(TLNeutrino + TLVisble).E()<<endl;
    
    for(int i = 0; i<vISR.size(); i++){
        MCParticle* a_MCP = vISR.at(i);
        v92.push_back(a_MCP);
    }
    
    std::vector<MCParticle*> vFL; vFL.clear();
    while( vFL.size() != numFL && num92 == 2 && count94 == 0){
        vFL.clear();
        for(int i = 0; i<v92.size(); i++){
            MCParticle* a_MCP = v92.at(i);
            if(a_MCP->getDaughters().size() != 0){
                for(int j = 0; j<a_MCP->getDaughters().size(); j++){
                    MCParticle* a_Dau = a_MCP->getDaughters()[j];
                    if(  abs(a_Dau->getPDG()) != 12 && abs(a_Dau->getPDG()) != 14 && abs(a_Dau->getPDG()) != 16 ){
                        vFL.push_back(a_Dau);
                    }
                }
            }
            else{
                vFL.push_back(a_MCP);
            }
        }
        
        TLorentzVector TLtest(0,0,0,0);
        ReconstructedParticleImpl* tempevent = new ReconstructedParticleImpl();
        for(int i = 0; i<vFL.size(); i++){
            MCParticle* a_MCP = vFL.at(i);
            TLorentzVector TLmcp(a_MCP->getMomentum(), a_MCP->getEnergy());
            TLtest += TLmcp;
            ReconstructedParticleImpl* b_Reco = new ReconstructedParticleImpl();
            b_Reco->setEnergy(a_MCP->getEnergy());
            b_Reco->setMomentum(a_MCP->getMomentum());
            tempevent->addParticle(b_Reco);
            LCRelationImpl *newrel = new LCRelationImpl(b_Reco, a_MCP, 1.0);
            otMCPSIMURelation->addElement( newrel );
        }

        otmcpsimufsp->addElement( tempevent );
        v92.clear();
        v92 = vFL;
        cout<<"vFL.size() : "<<vFL.size()<<endl;
        cout<<"TLtest.E() : "<<TLtest.E()<<endl;
    }


	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() << "   in run:  " << evt->getRunNumber() << std::endl ;
    evt->addCollection( otmcpsimufsp, _outmcpsimufsp.c_str() );
    evt->addCollection( otMCPSIMURelation, _outMCPSIMURelation.c_str() );
}

void MCParticleFSParticle::check( LCEvent * evt ) {
}

void MCParticleFSParticle::end() {

  std::cout << "MCParticleFSParticle::end()  " << std::endl;
}




std::vector<MCParticle*> getResult(std::vector<MCParticle*> vec, std::vector<MCParticle*> vec92){
    
    std::vector<MCParticle*> vecResult; vecResult.clear();
    
    for(int i = 0; i<vec.size(); i++){
        MCParticle* a_MCP = vec.at(i);
        if( std::find( vec92.begin(), vec92.end(), a_MCP ) != vec92.end() ){
            if( std::find(vecResult.begin(), vecResult.end(), a_MCP) == vecResult.end() ){
                vecResult.push_back(a_MCP);
            }
        }
        else{
            if(a_MCP->getParents().size() != 0){
                MCParticle* a_parent = a_MCP->getParents()[0];
                if( std::find( vec92.begin(), vec92.end(), a_parent ) != vec92.end() ){
                    if( std::find(vecResult.begin(), vecResult.end(), a_parent) == vecResult.end() ){
                        vecResult.push_back(a_parent);
                    }
                }
                else{
                    if( std::find(vecResult.begin(), vecResult.end(), a_parent) == vecResult.end() ){
                        vecResult.push_back(a_parent);
                    }
                }

            }
        }
    }
    
    return vecResult;
    
}


