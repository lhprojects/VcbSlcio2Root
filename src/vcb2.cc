#include <vcb2.hh>
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

vcb2 a_vcb2_instance;


vcb2::vcb2()
: Processor("vcb2"),
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
    
//    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
//                             "MCPSIMUFSP",
//                             "the final state particle after simulation",
//                             _outmcpsimufsp,
//                             std::string("ReconstructedParticle") );
//    
//    
//    registerOutputCollection( LCIO::LCRELATION,
//                             "mcpsimurelation",
//                             " relation between MCP and Reco after simulation",
//                             _outMCPSIMURelation,
//                             std::string("RelationusedSed") );
    
    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}


void findCombine(int index, vector<int> &a, vector<int>& tmp, const int& n, vector<vector<int> >& res){
	if(tmp.size() == n/2){
		res.push_back(tmp);
		return;
	}
	for(int i=index; i<=n; i++){
		tmp.push_back(a[i-1]);
		findCombine(i+1, a, tmp, n, res);
		tmp.pop_back();
	}
	return;
}
vector<vector<int> > Combine(int n, vector<int> a){
	vector<vector<int> > res;
	if(n<1 || n%2!=0 || a.size()<1)
		return res;
	vector<int> tmp;
	findCombine(1, a, tmp, n, res);
	return res;
}
vector<int> DifferenceSet(vector<int> total, vector<int> original, vector<int> left){
    for(int i = 0; i<total.size(); i++)
    {
        int element = total[i];
        vector<int >::iterator it;
        it = find(original.begin(), original.end(), element);
        if(it == original.end()){
            left.push_back(element);
        }
    }
    return left;
}
std::vector<vector<int> > pair4jets(int numjets);

static bool sortEn(ReconstructedParticle* a1, ReconstructedParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

void vcb2::init() {
    
    printParameters();
    TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));
    if (!tree_file->IsOpen()) {
        delete tree_file;
        tree_file=new TFile(_treeFileName.c_str(),"NEW");
    }
    
    _outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
    _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
    
    _outputTree->Branch("EventNr",       &eventNr,        "EventNr/I");
    _outputTree->Branch("Num",           &Num,            "Num/I");
    _outputTree->Branch("visEn",         &visEn,          "visEn/D");
    _outputTree->Branch("multiplicity",  &multiplicity,   "multiplicity/I");
    _outputTree->Branch("LeadElecEn",    &LeadElecEn,     "LeadElecEn/D");
    _outputTree->Branch("leadMuonEn",    &leadMuonEn,     "leadMuonEn/D");
    
    _outputTree->Branch("leadPionEn",    &leadPionEn,     "leadPionEn/D");
    _outputTree->Branch("leadGammaEn",    &leadGammaEn,     "leadGammaEn/D");


    _outputTree->Branch("subleadMuonEn", &subleadMuonEn,  "subleadMuonEn/D");
    _outputTree->Branch("EnA",           &EnA,            "EnA/D");

    _outputTree->Branch("leadMuonIMP",   &leadMuonIMP,    "leadMuonIMP/D");
    _outputTree->Branch("subleadMuonIMP",&subleadMuonIMP, "subleadMuonIMP/D");
    
    _outputTree->Branch("ratio15",       &ratio15,        "ratio15/D");
    _outputTree->Branch("ratio30",       &ratio30,        "ratio30/D");
    _outputTree->Branch("subratio15",    &subratio15,     "subratio15/D");
    _outputTree->Branch("subratio30",    &subratio30,     "subratio30/D");

    _outputTree->Branch("num_jet",       &num_jet,        "num_jet/I");
    _outputTree->Branch("JetsInvMass",   &JetsInvMass,    "JetsInvMass/D");
    _outputTree->Branch("JetsRecoilMass",&JetsRecoilMass, "JetsRecoilMass/D");
    _outputTree->Branch("jet1cosTheta",  &jet1cosTheta,   "jet1cosTheta/D");
    _outputTree->Branch("jet2cosTheta",  &jet2cosTheta,   "jet2cosTheta/D");
    _outputTree->Branch("jet1En",        &jet1En,         "jet1En/D");
    _outputTree->Branch("jet2En",        &jet2En,         "jet2En/D");
    _outputTree->Branch("jet1quark",     &jet1quark,      "jet1quark/I");
    _outputTree->Branch("jet2quark",     &jet2quark,      "jet2quark/I");

    _outputTree->Branch("leadMuonCharge",     &leadMuonCharge,      "leadMuonCharge/I");
    _outputTree->Branch("subleadMuonCharge",  &subleadMuonCharge,   "subleadMuonCharge/I");
    _outputTree->Branch("chargeA",            &chargeA,             "chargeA/I");


    _outputTree->Branch("num_quark",        &num_quark,        "num_quark/I");
    _outputTree->Branch("angle1",           &angle1,           "angle1/D");
    _outputTree->Branch("angle2",           &angle2,           "angle2/D");
    
    _outputTree->Branch("missPt",           &missPt,           "missPt/D");
    _outputTree->Branch("missM",            &missM,            "missM/D");
    _outputTree->Branch("angleLA",          &angleLA,          "angleLA/D");
    _outputTree->Branch("leadMcosTheta",    &leadMcosTheta,    "leadMcosTheta/D");

    
    
    _outputTree->Branch("tauDecay",          &tauDecay,            "tauDecay/I");

    _outputTree->Branch("HbbL",              &HbbL,                "HbbL[2]/F");
    _outputTree->Branch("HccL",              &HccL,                "HccL[2]/F");
    
    _outputTree->Branch("jet14m",              &jet14m,                "jet14m[4]/F");
    _outputTree->Branch("jet24m",              &jet24m,                "jet24m[4]/F");
    
    _outputTree->Branch("quark14m",              &quark14m,                "quark14m[4]/F");
    _outputTree->Branch("quark24m",              &quark24m,                "quark24m[4]/F");
    
    _outputTree->Branch("Y12",           &Y12,                 "Y12/D");
    _outputTree->Branch("Y23",           &Y23,                 "Y23/D");
    _outputTree->Branch("Y34",           &Y34,                 "Y34/D");
    
    _outputTree->Branch("quark1PDG",         &quark1PDG,           "quark1PDG/I");
    _outputTree->Branch("quark2PDG",         &quark2PDG,           "quark2PDG/I");
    _outputTree->Branch("quark3PDG",         &quark3PDG,           "quark3PDG/I");
    _outputTree->Branch("quark4PDG",         &quark4PDG,           "quark4PDG/I");

 
    Num = 0;
}

void vcb2::processEvent( LCEvent * evtP )
{
    
    if (evtP)
    {
        try{
            
            cout<<"Next Event *******************************************************************************************************"<<endl;
            eventNr = evtP->getEventNumber();
            cout<<"eventNr : "<<eventNr<<" Num : "<<Num<<endl;
            
            //ArborPFOs
            LCCollection* col_PFO = evtP->getCollection( "LICHPFOs" );
            int nPFO = col_PFO->getNumberOfElements();
            cout<<"nPFO : "<<nPFO<<endl;
            TLorentzVector TLPFO(0,0,0,0);
            
            
            std::vector<ReconstructedParticle*> vElec;      vElec.clear();
            std::vector<ReconstructedParticle*> vMuon;      vMuon.clear();
            std::vector<ReconstructedParticle*> vPion;      vPion.clear();
            std::vector<ReconstructedParticle*> vKaon;      vKaon.clear();
            std::vector<ReconstructedParticle*> vProton;      vProton.clear();
            std::vector<ReconstructedParticle*> vGamma;      vGamma.clear();



            for(int i = 0; i<nPFO; i++){
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                TLPFO += temp;
                int type = abs(pfo->getType());
                if(type == 11){vElec.push_back(pfo);}
                if(type == 13){vMuon.push_back(pfo);}
                if(type == 211){vPion.push_back(pfo);}
                if(type == 321){vKaon.push_back(pfo);}
                if(type == 2212){vProton.push_back(pfo);}
                if(type == 22){vGamma.push_back(pfo);}

            }
            
            cout<<"the total energy of ArborPFOs is : "<<TLPFO.E()<<endl;

            missPt = 999; missM = 999;
            TLorentzVector TLCM(0,0,0,240);
            missPt = ((TLCM - TLPFO).Vect()).Perp();
            missM = (TLCM - TLPFO).M();
            cout<<"missPt : "<<missPt<<endl;
            
            
           
                
            //the following code used to analysis the multiplicity information of final state particles
            visEn = 0; multiplicity = 0;
            visEn = TLPFO.E();
            multiplicity = nPFO;
            
            leadMuonEn = 0, LeadElecEn = 0, subleadMuonEn = 0, ratio15 = 0, ratio30 = 0, subratio15 = 0, subratio30 = 0;
            leadMuonCharge = 0, subleadMuonCharge = 0;
            leadMuonIMP = 0,    subleadMuonIMP = 0;
            leadMcosTheta = 0;

            
            sort(vElec.begin(), vElec.end(), sortEn);
            sort(vMuon.begin(), vMuon.end(), sortEn);
            sort(vPion.begin(), vPion.end(), sortEn);
            
            sort(vGamma.begin(), vGamma.end(), sortEn);

            

            ReconstructedParticle* leadMuon = NULL;
            ReconstructedParticle* leadMuonAround = NULL;
            ReconstructedParticle* subleadMuon = NULL;
            
            std::vector<ReconstructedParticle*> muonAround; muonAround.clear();
            
            

            ReconstructedParticle* leadElec = NULL;
            if(vMuon.size() > 0){
                leadMuon = vMuon.at(0);
                leadMuonEn = leadMuon->getEnergy();
                leadMuonCharge = leadMuon->getCharge();
                Track* leadMuonTrack = leadMuon->getTracks()[0];
                double pD0 = leadMuonTrack->getD0();
                double pZ0 = leadMuonTrack->getZ0();
                double psquare = pow( pD0*pD0 + pZ0*pZ0, 0.5 );
                leadMuonIMP = psquare;
                
                TLorentzVector TLseedMuon(leadMuon->getMomentum(), leadMuon->getEnergy());
                TVector3 TVseedMuon = TLseedMuon.Vect();
                leadMcosTheta = TVseedMuon.CosTheta();
                
                TLorentzVector TLaround15(0, 0, 0, 0);
                TLorentzVector TLaround30(0, 0, 0, 0);
                TLorentzVector subTLaround15(0, 0, 0, 0);
                TLorentzVector subTLaround30(0, 0, 0, 0);

                for(int i = 0; i<nPFO; i++){
                    ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                    TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                    TVector3 TVtemp = temp.Vect();
                    
                    if(pfo != leadMuon){
                        if( TVtemp.Angle(TVseedMuon) < 15./180 * 3.1415 ){ TLaround15 += temp; }
                        if( TVtemp.Angle(TVseedMuon) < 30./180 * 3.1415 ){ TLaround30 += temp;
                            if(abs(pfo->getType()) == 13){muonAround.push_back(pfo);}
                        }
                    }
                    
                }
                
                ratio15 = leadMuon->getEnergy()/( TLaround15.E() + leadMuon->getEnergy() );
                ratio30 = leadMuon->getEnergy()/( TLaround30.E() + leadMuon->getEnergy() );
                
                
                if(vMuon.size() > 1){
                    subleadMuon = vMuon.at(1);
                    subleadMuonEn = subleadMuon->getEnergy();
                    subleadMuonCharge = subleadMuon->getCharge();
                    
                    Track* subleadMuonTrack = subleadMuon->getTracks()[0];
                    double nD0 = subleadMuonTrack->getD0();
                    double nZ0 = subleadMuonTrack->getZ0();
                    double nsquare = pow( nD0*nD0 + nZ0*nZ0, 0.5 );
                    subleadMuonIMP = nsquare;
                    
                    for(int i = 0; i<nPFO; i++){
                        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                        TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                        TVector3 TVtemp = temp.Vect();
                        if(pfo != subleadMuon){
                            if( TVtemp.Angle(TVseedMuon) < 15./180 * 3.1415 ){ subTLaround15 += temp; }
                            if( TVtemp.Angle(TVseedMuon) < 30./180 * 3.1415 ){ subTLaround30 += temp; }
                        }
                    }
                    subratio15 = subleadMuon->getEnergy()/( subTLaround15.E() + subleadMuon->getEnergy() );
                    subratio30 = subleadMuon->getEnergy()/( subTLaround30.E() + subleadMuon->getEnergy() );
                }
            }
            
            
//            cout<<"vMuon.size() : "<<vMuon.size()<<" subleadMuonEn : "<<subleadMuonEn<<" subratio15 : "<<subratio15<<" subratio30 : "<<subratio30<<endl;
//            cout<<"leadMuonCharge : "<<leadMuonCharge<<" subleadMuonCharge : "<<subleadMuonCharge<<endl;
//            cout<<"leadMuonIMP : "<<leadMuonIMP<<" "<<subleadMuonIMP<<endl;
//            cout<<"missPt : "<<missPt<<" leadMcosTheta : "<<leadMcosTheta<<endl;
            
            if(vElec.size() > 0){
                leadElec = vElec.at(0);
                LeadElecEn = leadElec->getEnergy();
            }
            
            cout<<"vPion.size() : "<<vPion.size()<<endl;
            cout<<"vKaon.size() : "<<vKaon.size()<<endl;
            cout<<"vProton.size() : "<<vProton.size()<<endl;
            cout<<"vGamma.size() : "<<vGamma.size()<<endl;


            leadPionEn = 0;
            if(vPion.size() > 0){
                ReconstructedParticle* leadPion = vPion.at(0);
                leadPionEn = leadPion->getEnergy();
            }
            
            leadGammaEn = 0;
            if(vGamma.size() > 0){
                ReconstructedParticle* leadPion = vGamma.at(0);
                leadGammaEn = leadPion->getEnergy();
            }
            

            chargeA = 999; EnA = 999; angleLA = 999;
            if(muonAround.size() != 0){
                sort(muonAround.begin(), muonAround.end(), sortEn);
                leadMuonAround = muonAround.at(0);
                TLorentzVector TLLeadM(leadMuon->getMomentum(), leadMuon->getEnergy());
                TLorentzVector TLLA(leadMuonAround->getMomentum(), leadMuonAround->getEnergy());
                TVector3 TV3leadM = TLLeadM.Vect();
                TVector3 TV3LA = TLLA.Vect();
                angleLA = TV3leadM.Angle(TV3LA);
                EnA = leadMuonAround->getEnergy();
                chargeA = leadMuonAround->getCharge();
                cout<<"chargeA : "<<chargeA<<" EnA : "<<EnA<<" angleLA : "<<angleLA<<endl;
            }
            
            TLorentzVector TLT240(0, 0, 0, 240);
            
            tauDecay = 999;
            std::vector<MCParticle*> quarkvec; quarkvec.clear();
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            MCParticle* tauCandidate = NULL;
            
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents      = a_MCP->getParents().size();
                int NDaughters    = a_MCP->getDaughters().size();
                int PDG           = a_MCP->getPDG();
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TVector3 tempV3 = temp.Vect();
                
                if( abs(PDG) == 15){
                    
                    cout<<"find tau -------------------------------------------"<<endl;
                    
                    for(int j = 0; j<a_MCP->getDaughters().size(); j++){
                        MCParticle* tauDau = a_MCP->getDaughters()[j];
                        if(abs(tauDau->getPDG()) == 16){
                            tauCandidate = a_MCP;
                            break;
                        }
                    }
                }
                
                if(NParents == 0 && (abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6) ){
                    quarkvec.push_back(a_MCP);
                    
                }
                
                if( PDG == 25 && NDaughters == 4 ){
                    cout<<"quarkvec.size() : "<<quarkvec.size()<<endl;
                    
                    MCParticle* HDau = a_MCP->getDaughters()[0];
                    int HDauPDG = HDau->getPDG();
                    if( abs(HDauPDG) == 1 || abs(HDauPDG) == 2 || abs(HDauPDG) == 3 || abs(HDauPDG) == 4 || abs(HDauPDG) == 5 || abs(HDauPDG) == 6 || abs(HDauPDG) == 21 ){
                        
                        cout<<"HDauPDG : "<<HDauPDG<<endl;
                        
                        quarkvec.push_back( HDau );
                        MCParticle* HDau2 = a_MCP->getDaughters()[1];
                        quarkvec.push_back( HDau2 );
                        cout<<"quarkvec.size() : "<<quarkvec.size()<<endl;
                    }
                }
            }
            
            if(tauCandidate != NULL){
                for(int i = 0; i<tauCandidate->getDaughters().size(); i++){
                    MCParticle* tauD = tauCandidate->getDaughters()[i];
                    if( abs(tauD->getPDG()) == 11 || abs(tauD->getPDG()) == 13 ){ tauDecay = tauD->getPDG(); }
                }
                if( tauDecay == 999  ){ tauDecay = 1; }
            }
            
            cout<<"tauDecay : "<<tauDecay<<endl;
            
            for(int i = 0; i<2; i++){
                HbbL[i] = 999;
                HccL[i] = 999;
            }
            
            for(int i = 0; i<4; i++){
                jet14m[i] = 999;
                jet24m[i] = 999;
                quark14m[i] = 999;
                quark24m[i] = 999;
            }

            LCCollection* col_Jet = evtP->getCollection( "RefinedJets" );
            num_jet = col_Jet->getNumberOfElements();
            cout<<"num_jet : "<<num_jet<<endl;
            JetsInvMass = 0;
            jet1cosTheta = 999, jet2cosTheta = 999, jet1En = 999, jet2En = 999;
            TLorentzVector TLJets(0,0,0,0);
            
            Y12 = 999, Y23 = 999, Y34 = 999;

            
            if(num_jet >= 1){

                PIDHandler pidh(col_Jet);
                int algo   = 0;
                int ibtag  =-1;
                int ictag  =-1;
                int ibctag =-1;
                int icat   =-1;
                
                
                
                algo = pidh.getAlgorithmID("yth");
                int iy12 = -1, iy23 = -1, iy34 = -1;
                iy12 = pidh.getParameterIndex (algo, "y12" );
                iy23 = pidh.getParameterIndex (algo, "y23" );
                iy34 = pidh.getParameterIndex (algo, "y34" );
                ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
                const ParticleID &pid = pidh.getParticleID(jet, algo);
                Y12 = pid.getParameters()[iy12];
                Y23 = pid.getParameters()[iy23];
                Y34 = pid.getParameters()[iy34];

                
                
                algo      = pidh.getAlgorithmID("lcfiplus");
                ibtag     = pidh.getParameterIndex (algo,  "BTag");
                ictag     = pidh.getParameterIndex (algo,  "CTag");
                ibctag    = pidh.getParameterIndex (algo,  "BCTag");
                icat      = pidh.getParameterIndex (algo,  "Category");
                
                
                for(int i = 0; i<num_jet; i++){
                    ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(i));
                    TLorentzVector TLtemp(jet->getMomentum(), jet->getEnergy());
                    TLJets += TLtemp;
                    const ParticleID &pid = pidh.getParticleID(jet, algo);
                    double btag  = pid.getParameters()[ibtag];
                    double ctag  = pid.getParameters()[ictag];
                    double bctag = pid.getParameters()[ibctag];
                    double cat   = pid.getParameters()[icat];
                    double otag = 1 - btag - ctag;
                    
                }
            }
            
            
            JetsInvMass = TLJets.M();
            JetsRecoilMass = 0;
            JetsRecoilMass = (TLT240 - TLJets).M();
            
            jet1quark = 0, jet2quark = 0;
            num_quark = 0;
            num_quark = quarkvec.size();
            angle1 = 999, angle2 = 999;
            
            quark2PDG = 0, quark3PDG = 0, quark4PDG = 0, quark1PDG = 0;
            
            if(num_jet == 2){
                
                PIDHandler pidh(col_Jet);
                int algo   = 0;
                int ibtag  =-1;
                int ictag  =-1;
                int ibctag =-1;
                int icat   =-1;
                algo      = pidh.getAlgorithmID("lcfiplus");
                ibtag     = pidh.getParameterIndex (algo,  "BTag");
                ictag     = pidh.getParameterIndex (algo,  "CTag");
                ibctag    = pidh.getParameterIndex (algo,  "BCTag");
                icat      = pidh.getParameterIndex (algo,  "Category");
                
                ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
                
                const ParticleID &pid1 = pidh.getParticleID(jet1, algo);
                double btag1  = pid1.getParameters()[ibtag];
                double ctag1  = pid1.getParameters()[ictag];

                
                TLorentzVector TLjet1(jet1->getMomentum(), jet1->getEnergy());
                TVector3 TVjet1 = TLjet1.Vect();
                jet1cosTheta = TVjet1.CosTheta();
                jet1En = jet1->getEnergy();
                jet14m[0] = jet1->getMomentum()[0]; jet14m[1] = jet1->getMomentum()[1];
                jet14m[2] = jet1->getMomentum()[2]; jet14m[3] = jet1->getEnergy();

                ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(1));
                
                const ParticleID &pid2 = pidh.getParticleID(jet2, algo);
                double btag2  = pid2.getParameters()[ibtag];
                double ctag2  = pid2.getParameters()[ictag];

                TLorentzVector TLjet2(jet2->getMomentum(), jet2->getEnergy());
                TVector3 TVjet2 = TLjet2.Vect();
                jet2cosTheta = TVjet2.CosTheta();
                jet2En = jet2->getEnergy();
                
                jet24m[0] = jet1->getMomentum()[0]; jet24m[1] = jet1->getMomentum()[1];
                jet24m[2] = jet1->getMomentum()[2]; jet24m[3] = jet1->getEnergy();
                
                cout<<"btag1 : "<<btag1<<" btag2 "<<btag2<<endl;
                cout<<"jet1cosTheta : "<<jet1cosTheta<<" "<<jet1En<<" "<<jet2cosTheta<<" "<<jet2En<<endl;
                
                HbbL[0] = btag1; HbbL[1] = btag2;
                HccL[0] = ctag1; HccL[1] = ctag2;
                
                if(quarkvec.size() > 1){
                    MCParticle* quark1 = quarkvec.at(0);
                    MCParticle* quark2 = quarkvec.at(1);
                    TLorentzVector TLquark1(quark1->getMomentum(), quark1->getEnergy());
                    TLorentzVector TLquark2(quark2->getMomentum(), quark2->getEnergy());
                    TVector3 TVquark1 = TLquark1.Vect();
                    TVector3 TVquark2 = TLquark2.Vect();
                
                    if(TVjet1.Angle(TVquark1) + TVjet2.Angle(TVquark2) < TVjet1.Angle(TVquark2) + TVjet2.Angle(TVquark1)){
                        jet1quark = quark1->getPDG();
                        jet2quark = quark2->getPDG();
                        angle1 = TVjet1.Angle(TVquark1);
                        angle2 = TVjet2.Angle(TVquark2);
                        
                        quark14m[0] = quark1->getMomentum()[0]; quark14m[1] = quark1->getMomentum()[1];
                        quark14m[2] = quark1->getMomentum()[2]; quark14m[3] = quark1->getEnergy();
                        quark24m[0] = quark2->getMomentum()[0]; quark24m[1] = quark2->getMomentum()[1];
                        quark24m[2] = quark2->getMomentum()[2]; quark24m[3] = quark2->getEnergy();
                    }
                    else{
                        jet1quark = quark2->getPDG();
                        jet2quark = quark1->getPDG();
                        angle1 = TVjet1.Angle(TVquark2);
                        angle2 = TVjet2.Angle(TVquark1);
                        
                        quark14m[0] = quark2->getMomentum()[0]; quark14m[1] = quark2->getMomentum()[1];
                        quark14m[2] = quark2->getMomentum()[2]; quark14m[3] = quark2->getEnergy();
                        quark24m[0] = quark1->getMomentum()[0]; quark24m[1] = quark1->getMomentum()[1];
                        quark24m[2] = quark1->getMomentum()[2]; quark24m[3] = quark1->getEnergy();
                    }
                
                    quark1PDG = quark1->getPDG();
                    quark2PDG = quark2->getPDG();
                    
                    if(quarkvec.size() > 2){
                        MCParticle* quark3 = quarkvec.at(2);
                        MCParticle* quark4 = quarkvec.at(3);
                        quark3PDG = quark3->getPDG();
                        quark4PDG = quark4->getPDG();
                    }
                }
            }
            

        }catch (lcio::DataNotAvailableException err) {  }

    }
    
    _outputTree->Fill();
    Num ++;
}



void vcb2::end()
{
    
    if (_outputTree) {
        
        TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        //tree_file->cd();
        tree_file->Write();
        delete tree_file;
        //tree_file->Close();
    }
    
}





