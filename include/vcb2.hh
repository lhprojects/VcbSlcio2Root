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
    
    
        double ratio15, ratio30, subratio15, subratio30;
        double visEn, missPt, missM;
        double leadMuonEn,  LeadElecEn, subleadMuonEn, leadPionEn, leadGammaEn;
        int    multiplicity;
    
    int tauDecay;
    
    int jet1quark, jet2quark, num_jet, num_quark;
    double angle1, angle2;
    
    double JetsInvMass, JetsRecoilMass;
    float HbbL[2], HccL[2];
    
    float jet14m[4], jet24m[4];
    float quark14m[4], quark24m[4];
    double Y12, Y23, Y34;

    
    double jet1cosTheta, jet2cosTheta, jet1En, jet2En;
    
    int subleadMuonCharge, leadMuonCharge, chargeA;
    
    double leadMuonIMP, subleadMuonIMP, leadMcosTheta;
    
    int quark1PDG, quark2PDG, quark3PDG, quark4PDG;
    double EnA, angleLA;
    

		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


