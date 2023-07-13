#ifndef _MC2Reconstructed_hh_
#define _MC2Reconstructed_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

class TTree;

class MC2Reconstructed  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new MC2Reconstructed ; }

		MC2Reconstructed();

		~MC2Reconstructed() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:

		int _overwrite;
    
        std::string _outmcpsimufsp;
        std::string _outMCPSIMURelation;
    
};

#endif


