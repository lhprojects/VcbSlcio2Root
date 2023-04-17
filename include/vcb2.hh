#ifndef _vcb2_hh_
#define _vcb2_hh_

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

class vcb2  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new vcb2 ; }

		vcb2();

		~vcb2() {};

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
    
    float HbbL[2], HccL[2];
    
    float jet14m[4],   jet24m[4];
    float quark14m[4], quark24m[4];
    float lead1HadEn, lead2HadEn;
    int   lead1PID,   lead2PID;
    float lead1EndP[3], lead2EndP[3];
    int num_92BDaus;
    
    double thrust, Y12, Y23, Y34;
    
    int quark1PDG, quark2PDG;

    

		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif

