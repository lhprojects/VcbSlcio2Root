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
		void setBranchAndValue(TTree *_outputTree_);

		void processEvent( LCEvent * evtP );
		void doProcessEvent( LCEvent * evtP );
		void fillJets(LCCollection *col_Jet, std::vector<MCParticle*> &quarkvec);
		void fillLeptons(std::vector<MCParticle*> &leptonvec);
		ReconstructedParticle *saveIsoLepton(std::vector<ReconstructedParticle *> const &isoLeps, LCCollection *col_PFO);
		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
//		TFile *tree_file;

		int _overwrite;
		TTree *_outputTree;

		int eventNr;
		int Num;
		double _centerOfMassEnergy;

		int _isoLepPDG;
    
        double visEn, missPt, missM;
        int    multiplicity;

        double leadMuonEn,  leadElecEn, leadPionEn, leadGammaEn;

        double leadLepEn;
        double leadLepCharge;
		double leadD0;
		double leadZ0;
		double leadIMP;
		double leadLepRatio15;
		double leadLepRatio30;
		double leadLepCostheta;

        double subleadLepEn;
        double subleadLepCharge;
		double subleadD0;
		double subleadZ0;
		double subleadIMP;
		double subleadLepRatio15;
		double subleadLepRatio30;
		double subleadLepCostheta;

        double leadConeSubEn;
		double leadConeAngle;
		double leadConeTotalMass;
		double leadConeTotalRatio;
	

		double sametypePairMass;
		double difftypePairMass;
		double sametypePairRecoilMass;
		double difftypePairRecoilMass;


    int tauDecay;
    
    int jet1quark, jet2quark, num_jet;
    
    double JetsInvMass, JetsRecoilMass;
    double jet1cosTheta, jet2cosTheta, jet1En, jet2En;
    float HbbL[2], HccL[2];
	int jet1MCPDG, jet2MCPDG;
    float jet14m[4], jet24m[4];
    double Y12, Y23, Y34;

	int num_quark;
    double angle1, angle2;
    float quark14m[4], quark24m[4];

	int num_lepton;
    float lepton14m[4], lepton24m[4];
	int lepton1PDG, lepton2PDG;
    
    
    int subleadMuonCharge, leadMuonCharge, chargeA;
    
    double leadMuonIMP, subleadMuonIMP, leadMcosTheta;
    
    int quark1PDG, quark2PDG, quark3PDG, quark4PDG;    

	std::string _fileName;
	std::ostream *_output;
	std::string _histFileName;
};

#endif


