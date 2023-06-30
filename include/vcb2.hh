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
#include "TLorentzVector.h"

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
		void fracInCone(ReconstructedParticle *part, LCCollection *col_PFO, bool leading);
		void fillParticles(LCCollection *col_PFO,
			std::vector<ReconstructedParticle *> &vElec,
			std::vector<ReconstructedParticle *> &vMuon,
			std::vector<ReconstructedParticle *> &vPion,
			std::vector<ReconstructedParticle *> &vKaon,
			std::vector<ReconstructedParticle *> &vProton,
			std::vector<ReconstructedParticle *> &vGamma,
			TLorentzVector &TLPFO
			);

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
    
        float visEn, missPt, missM;
        int    multiplicity;

		int nPFOs;
		int nGoodPFOs;
		int nChargedPFOs;
		int nGoodChargedPFOs;

        float leadMuonEn,  leadElecEn, leadPionEn, leadGammaEn;

        float leadLepEn;
        float leadLepCharge;
		float leadD0;
		float leadZ0;
		float leadIMP;
		float leadLepRatio15;
		float leadLepRatio30;

		float leadLepEn_501[5];
		float leadLepEn_22[5];
		float leadLepEn_21120[5];
		float leadLepEn_other[5];

		float subleadLepEn_501[5];
		float subleadLepEn_22[5];
		float subleadLepEn_21120[5];
		float subleadLepEn_other[5];

		float leadLepCostheta;

        float subleadLepEn;
        float subleadLepCharge;
		float subleadD0;
		float subleadZ0;
		float subleadIMP;
		float subleadLepRatio15;
		float subleadLepRatio30;
		float subleadLepCostheta;

		int leadConePDG;
        float leadConeSubEn;
		float leadConeAngle;
		float leadConeTotalMass;
		float leadConeTotalRatio;
	

		float sametypePairMass;
		float difftypePairMass;
	    float sametypePairRecoilMass;
		float difftypePairRecoilMass;



    int tauDecay;
    
    int jet1quark, jet2quark, num_jet;
    
    float JetsInvMass, JetsRecoilMass;
    float jet1cosTheta, jet2cosTheta, jet1En, jet2En;
    float HbbL[2], HccL[2];
	int jet1MCPDG, jet2MCPDG;
    float jet14m[4], jet24m[4];
    float Y12, Y23, Y34;

	float leadLepM4[4];
	float subleadLepM4[4];

	int num_quark;
    float angle1, angle2;
    float quark14m[4], quark24m[4];

	int num_lepton;
    float lepton14m[4], lepton24m[4];
	int lepton1PDG, lepton2PDG;
    
    
    int subleadMuonCharge, leadMuonCharge, chargeA;
    
    float leadMuonIMP, subleadMuonIMP, leadMcosTheta;
    
    int quark1PDG, quark2PDG, quark3PDG, quark4PDG;    

	std::string _fileName;
	std::ostream *_output;
	std::string _histFileName;
};

#endif


