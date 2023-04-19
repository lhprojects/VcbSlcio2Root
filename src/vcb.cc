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

    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "MCPSIMUFSP2",
                            "the final state particle after simulation",
                            _outmcpsimufsp2,
                            std::string("ReconstructedParticle") );
    
    registerOutputCollection( LCIO::LCRELATION,
                             "mcpsimurelation",
                             " relation between MCP and Reco after simulation",
                             _outMCPSIMURelation,
                             std::string("RelationusedSed") );
    
        registerOutputCollection( LCIO::LCRELATION,
                             "mcpsimurelation2",
                             " relation between MCP and Reco after simulation",
                             _outMCPSIMURelation2,
                             std::string("RelationusedSed") );

    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}




// from large to small
struct SortEn {
    bool operator()(ReconstructedParticle const * a1, ReconstructedParticle const * a2) const{
        return a1->getEnergy() > a2->getEnergy();
    }
};

static SortEn sortEn;

void vcb::init() {    
    printParameters();
    Num = 0;
}


ReconstructedParticle *clone(ReconstructedParticle * pfo)
{
    if(false) {
        // yong feng's code seems also working
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
            for(unsigned j = 0; j<pfo->getClusters().size(); j++){
                a_Reco->addCluster( pfo->getClusters()[j] );
            }
        }
        return a_Reco;
    } else if(false) {
        ReconstructedParticleImpl* a_Reco = new ReconstructedParticleImpl();
        // fields: need to copy
        // int _type{0} ;
        // double _momentum[3] = {0.,0.,0.} ;
        // double _energy{0.} ;
        // EVENT::FloatVec _cov{} ;
        // double _mass{0.} ;
        // float _charge{} ;
        // float _reference[3] = {0.,0.,0.}  ;
        // EVENT::ParticleID* _pidUsed{ NULL} ;
        // float _goodnessOfPID{0.} ;
        // EVENT::ParticleIDVec _pid{} ;
        // EVENT::ReconstructedParticleVec _particles{} ;
        // EVENT::ClusterVec _clusters{} ;
        // EVENT::TrackVec _tracks{} ;
        // EVENT::Vertex* _sv{} ;


        a_Reco->setType( pfo->getType() );
        a_Reco->setMomentum( pfo->getMomentum() );
        a_Reco->setEnergy( pfo->getEnergy() );
        a_Reco->setCovMatrix(pfo->getCovMatrix());
        a_Reco->setMass(pfo->getMass());
        a_Reco->setCharge(pfo->getCharge());
        a_Reco->setReferencePoint (pfo->getReferencePoint());
        a_Reco->setParticleIDUsed(pfo->getParticleIDUsed());
        a_Reco->setGoodnessOfPID(pfo->getGoodnessOfPID());
        
        for(unsigned i = 0; i < pfo->getParticleIDs().size(); ++i) {
            a_Reco->addParticleID(pfo->getParticleIDs()[i]);
        }

        for(unsigned i = 0; i < pfo->getParticles().size(); ++i) {
            a_Reco->addParticle(pfo->getParticles()[i]);
        }

        for(unsigned i = 0; i < pfo->getClusters().size(); ++i) {
            a_Reco->addCluster(pfo->getClusters()[i]);
        }

        for(unsigned i = 0; i < pfo->getTracks().size(); ++i) {
            a_Reco->addTrack(pfo->getTracks()[i]);
        }
        a_Reco->setStartVertex(a_Reco->getStartVertex());

        return a_Reco;
    } else{
        // do we really need to copy the obj?
        return pfo;
    }
}

struct Copy {

    LCCollectionVec* otmcpsimufsp;
    LCCollectionVec* otMCPSIMURelation;
    LCFlagImpl newMCPSIMUlinkflag;
    std::string _outmcpsimufsp;
    std::string _outMCPSIMURelation;
    int NumPart;

    Copy(std::string _outmcpsimufsp, std::string _outMCPSIMURelation) :
        _outmcpsimufsp(_outmcpsimufsp),
        _outMCPSIMURelation(_outMCPSIMURelation) {

        otmcpsimufsp = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
        otMCPSIMURelation = new LCCollectionVec( LCIO::LCRELATION );
        newMCPSIMUlinkflag.setBit(LCIO::CHBIT_LONG);
        otMCPSIMURelation->setFlag(newMCPSIMUlinkflag.getFlag());
        NumPart = 0;
    }

    void addPart(ReconstructedParticle *pfo) {
        //pfo = clone(pfo);
        ReconstructedParticle *a_Reco = clone(pfo);
        otmcpsimufsp->addElement(a_Reco );
        LCRelationImpl *newrel = new LCRelationImpl(a_Reco, pfo, 1.0);
        otMCPSIMURelation->addElement( newrel );
        NumPart += 1;
    }

    void added2Evt(LCEvent * evtP) {

        printf("%d parts added to %s\n", NumPart, _outmcpsimufsp.c_str());
        evtP->addCollection(otmcpsimufsp, _outmcpsimufsp.c_str() );
        evtP->addCollection(otMCPSIMURelation, _outMCPSIMURelation.c_str() );
    }

};

void vcb::processEvent( LCEvent * evtP )
{
    
    if (evtP)
    {
        Copy exclMuon(_outmcpsimufsp, _outMCPSIMURelation);
        Copy exclElec(_outmcpsimufsp2, _outMCPSIMURelation2);

        try {
            
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int nPFO = col_PFO->getNumberOfElements();
            cout<<"nPFO : "<< nPFO <<endl;
            
            
            std::vector<ReconstructedParticle*> vElec;
            std::vector<ReconstructedParticle*> vMuon;

            for(int i = 0; i < nPFO; i++) {
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                
                int type = abs(pfo->getType());
                if(type == 11){ vElec.push_back(pfo); }
                if(type == 13){ vMuon.push_back(pfo); }
                else {
                    exclMuon.addPart(pfo);
                    exclElec.addPart(pfo);
                } 
            }
            

            std::sort(vElec.begin(), vElec.end(), sortEn);
            std::sort(vMuon.begin(), vMuon.end(), sortEn);
            
            for (int i = 0; i < (int)vMuon.size(); i++)
            {
                ReconstructedParticle *pfo = vMuon.at(i);
                if(i != 0) {
                    exclMuon.addPart(pfo);
                }
                exclElec.addPart(pfo);                
            }

            for (int i = 0; i < (int)vElec.size(); i++)
            {
                ReconstructedParticle *pfo = vElec.at(i);
                if(i != 0) {
                    exclElec.addPart(pfo);
                }
                exclMuon.addPart(pfo);  
            }           

            exclMuon.added2Evt(evtP);
            exclElec.added2Evt(evtP);

        } catch (lcio::DataNotAvailableException err) {
            printf("process %dth event error", Num);
        }
    }
    
    Num++;
}

void vcb::end()
{    
}


