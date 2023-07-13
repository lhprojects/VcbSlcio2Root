#include <MC2Reconstructed.hh>
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


using namespace std;

MC2Reconstructed a_MC2Reconstructed_instance;


MC2Reconstructed::MC2Reconstructed()
: Processor("MC2Reconstructed")
{
    _description = "MC2Reconstructed" ;

    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             "ReconstructedParticleCollectionName",
                             "the final state particle after simulation",
                             _outmcpsimufsp,
                             std::string("ReconstructedParticle") );
    
    
    registerOutputCollection( LCIO::LCRELATION,
                             "ReconstructedParticleRelationCollectionName",
                             " relation between MCP and Reco after simulation",
                             _outMCPSIMURelation,
                             std::string("RelationusedSed") );
    
    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}

void MC2Reconstructed::init() {    
    printParameters();
}

void MC2Reconstructed::processEvent( LCEvent * evtP )
{
        
    LCCollectionVec* otmcpsimufsp = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
    LCCollectionVec* otMCPSIMURelation = new LCCollectionVec( LCIO::LCRELATION );
    LCFlagImpl newMCPSIMUlinkflag;
    newMCPSIMUlinkflag.setBit(LCIO::CHBIT_LONG);
    otMCPSIMURelation->setFlag(newMCPSIMUlinkflag.getFlag());
    
    if (evtP)
    {
        try{
            
            LCCollection* col_PFO = evtP->getCollection( "MCParticle" );
            int nPFO = col_PFO->getNumberOfElements();
            
            
            for(int i = 0; i<nPFO; i++){
                MCParticle* pfo = dynamic_cast<MCParticle*>(col_PFO->getElementAt(i));

                int type = abs(pfo->getPDG());
                if(pfo->getGeneratorStatus() == 1 && type != 12 && type != 14 && type != 16)
                {
                    if(type == 22) {
                        if(pfo->getEnergy() < 1E-3){
                            continue;
                        }
                    }
                    
                    TVector3 p3(pfo->getMomentum());
                    double costheta = p3.CosTheta();
                    if(fabs(costheta) > 0.995) {
                        continue;
                    }

                    ReconstructedParticleImpl* a_Reco = new ReconstructedParticleImpl();
                    a_Reco->setEnergy( pfo->getEnergy() );
                    a_Reco->setMomentum( pfo->getMomentum() );
                    a_Reco->setType( pfo->getPDG() );                    
                    a_Reco->setCharge(pfo->getCharge());
                    otmcpsimufsp->addElement( a_Reco );
                    LCRelationImpl *newrel = new LCRelationImpl(a_Reco, pfo, 1.0);
                    otMCPSIMURelation->addElement( newrel );
                }                
            }
            //printf("nphoton %d\n", nphoton);

            
        }catch (lcio::DataNotAvailableException err) {  }

    }
    
    evtP->addCollection( otmcpsimufsp, _outmcpsimufsp.c_str() );
    evtP->addCollection( otMCPSIMURelation, _outMCPSIMURelation.c_str() );
}



void MC2Reconstructed::end()
{
}


