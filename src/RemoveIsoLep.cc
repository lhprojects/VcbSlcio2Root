#include <RemoveIsoLep.hh>
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

RemoveIsoLep a_vcb_instance;


RemoveIsoLep::RemoveIsoLep()
: Processor("RemoveIsoLep"),
_output(0)
{
    _description = "Print MC Truth" ;
    
    
    
    _isoLepPDG = 13;
    registerProcessorParameter("IsoLepPDG",
                               "IsoLepPDG",
                               _isoLepPDG,
                               _isoLepPDG);

     registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             "InputCollectionName",
                             "input the final state particle after simulation",
                             _inputCollectionName,
                             std::string("ArborFPOs") );

   
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



// from large to small
struct SortEn {
    bool operator()(ReconstructedParticle const * a1, ReconstructedParticle const * a2) const{
        if(!a1) cout << "a1 NULL" << endl;
        if(!a2) cout << "a2 NULL" << endl;
        return a1->getEnergy() > a2->getEnergy();
    }
};

static SortEn sortEn;

void RemoveIsoLep::init() {    
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
    } else if(true) {
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
        
        // cause problem sometimes ?
        // need deep copy
        //for(unsigned i = 0; i < pfo->getParticleIDs().size(); ++i) {
        //    a_Reco->addParticleID(pfo->getParticleIDs()[i]);
        //}

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

void RemoveIsoLep::processEvent( LCEvent * evtP )
{
    
    if (evtP)
    {
        Copy exclMuon(_outmcpsimufsp, _outMCPSIMURelation);
        //Copy exclElec(_outmcpsimufsp2, _outMCPSIMURelation2);

        try {
            
            LCCollection* col_PFO = evtP->getCollection( _inputCollectionName );
            int nPFO = col_PFO->getNumberOfElements();
            cout<<"nPFO : "<< nPFO <<endl;
            

            ReconstructedParticle* leadIsoLep = NULL;
            double leadIsoLepEn = 0;
            for(int i = 0; i < nPFO; i++) {
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                
                if(_isoLepPDG == abs(pfo->getType()) ) {
                    if(pfo->getEnergy() > leadIsoLepEn) {
                        leadIsoLepEn = pfo->getEnergy();
                        leadIsoLep = pfo;
                        //std::cout << leadIsoLepEn << std::endl;
                    }
                }
            }
    
            for(int i = 0; i < nPFO; i++) {
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                if(pfo != leadIsoLep) {
                    exclMuon.addPart(pfo);                    
                }
            }

            exclMuon.added2Evt(evtP);
            //exclElec.added2Evt(evtP);

        } catch (lcio::DataNotAvailableException err) {
            printf("process %dth event error", Num);
        }
    }
    
    Num++;
}

void RemoveIsoLep::end()
{    
}


