#ifndef _vcb_hh_
#define _vcb_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

class TTree;

class vcb  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new vcb ; }

		vcb();

		~vcb() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
//		TFile *tree_file;

		int _overwrite;
		TTree *_outputTree;

		unsigned int eventNr;
		int Num;
//    
//    
//        double ratio;
//        double visEn;
//        double leadMuonEn,  LeadElecEn, subleadMuonEn;
//
//        int    multiplicity;

    
        std::string _outmcpsimufsp;
        std::string _outMCPSIMURelation;
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


