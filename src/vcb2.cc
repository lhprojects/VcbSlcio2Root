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
    _description = "Print MC Truth";

    _treeFileName = "MCTruth.root";
    registerProcessorParameter("TreeOutputFile",
                               "The name of the file to which the ROOT tree will be written",
                               _treeFileName,
                               _treeFileName);

    _treeName = "Evts";
    registerProcessorParameter("TreeName",
                               "The name of the ROOT tree",
                               _treeName,
                               _treeName);

    _isoLepPDG = 13;
    registerProcessorParameter("IsoLepPDG",
                               "IsoLepPDG",
                               _isoLepPDG,
                               _isoLepPDG);

    _centerOfMassEnergy = 240;
    registerProcessorParameter("CenterOfMassEnergy",
                               "The name of the ROOT tree",
                               _centerOfMassEnergy,
                               _centerOfMassEnergy);

    _overwrite = 1;
    registerProcessorParameter("OverwriteFile",
                               "If zero an already existing file will not be overwritten.",
                               _overwrite,
                               _overwrite);
}

TLorentzVector lorentzV4(ReconstructedParticle *mcp) {
    return TLorentzVector(mcp->getMomentum(), mcp->getEnergy());
}

TLorentzVector lorentzV4(MCParticle *mcp) {
    return TLorentzVector(mcp->getMomentum(), mcp->getEnergy());
}


static bool sortEn(ReconstructedParticle *a1, ReconstructedParticle *a2)
{
    return a1->getEnergy() >= a2->getEnergy();
}

struct TT
{
    TT(TTree *_outputTree) : _outputTree(_outputTree) { }
    TTree *_outputTree;

    template<class T>
    void Branch(char const *name, T const &value, char const *id) {
        if(_outputTree) {
            _outputTree->Branch(name, value, id);
        }
    }

};

void vcb2::setBranchAndValue(TTree *_outputTree_) {

    TT att(_outputTree_);
    TT *_outputTree = &att;

    eventNr =-1;
    // don't set Num
    _outputTree->Branch("EventNr", &eventNr, "EventNr/I"); 
    _outputTree->Branch("Num", &Num, "Num/I");

    multiplicity = 0;
    nPFOs = 0;
    nGoodPFOs = 0;
    nChargedPFOs = 0;
    nGoodChargedPFOs = 0;

    visEn = 0;
    missPt = 0;
    missM = 0;

    _outputTree->Branch("multiplicity", &multiplicity, "multiplicity/I");
    _outputTree->Branch("nPFOs", &nPFOs, "nPFOs/I");
    _outputTree->Branch("nGoodPFOs", &nGoodPFOs, "nGoodPFOs/I");
    _outputTree->Branch("nChargedPFOs", &nChargedPFOs, "nChargedPFOs/I");
    _outputTree->Branch("nGoodChargedPFOs", &nGoodChargedPFOs, "nGoodChargedPFOs/I");

    _outputTree->Branch("visEn", &visEn, "visEn/F");
    _outputTree->Branch("missPt", &missPt, "missPt/F");
    _outputTree->Branch("missM", &missM, "missM/F");


    _outputTree->Branch("LeadElecEn", &leadElecEn, "LeadElecEn/F"); leadElecEn = 0;
    _outputTree->Branch("leadMuonEn", &leadMuonEn, "leadMuonEn/F"); leadMuonEn = 0;
    _outputTree->Branch("leadPionEn", &leadPionEn, "leadPionEn/F"); leadPionEn = 0;
    _outputTree->Branch("leadGammaEn", &leadGammaEn, "leadGammaEn/F"); leadGammaEn = 0;


    for(int i = 0; i < 4; ++i) leadLepM4[i] = 0;
    _outputTree->Branch("leadLepM4", &leadLepM4, "leadLepM4[4]/F");
    _outputTree->Branch("leadLepEn", &leadLepEn, "leadLepEn/F"); leadLepEn = 0;
    _outputTree->Branch("leadLepCharge", &leadLepCharge, "leadLepCharge/F"); leadLepCharge = 999;
    _outputTree->Branch("leadD0", &leadD0, "leadD0/F"); leadD0 = 0;
    _outputTree->Branch("leadZ0", &leadZ0, "leadZ0/F"); leadZ0 = 0;
    _outputTree->Branch("leadIMP", &leadIMP, "leadIMP/F"); leadIMP = 0;
    _outputTree->Branch("leadLepCostheta", &leadLepCostheta, "leadLepCostheta/F"); leadLepCostheta= 999;
    _outputTree->Branch("leadLepRatio15", &leadLepRatio15, "leadLepRatio15/F"); leadLepRatio15 = 0;
    _outputTree->Branch("leadLepRatio30", &leadLepRatio30, "leadLepRatio30/F"); leadLepRatio30 = 0;

    for(int i = 0; i < 4; ++i) subleadLepM4[i] = 0;

    _outputTree->Branch("subleadLepM4", &subleadLepM4, "subleadLepM4[4]/F");

    _outputTree->Branch("subleadLepEn", &subleadLepEn, "subleadLepEn/F"); subleadLepEn = 0;
    _outputTree->Branch("subleadLepCharge", &subleadLepCharge, "subeadLepCharge/F"); subleadLepCharge = 0;
    _outputTree->Branch("subleadD0", &subleadD0, "subleadD0/F"); subleadD0 = 0;
    _outputTree->Branch("subleadZ0", &subleadZ0, "subleadZ0/F"); subleadZ0 = 0;
    _outputTree->Branch("subleadIMP", &subleadIMP, "subleadIMP/F"); subleadIMP = 0;
    _outputTree->Branch("subleadLepCostheta", &subleadLepCostheta, "subleadLepCostheta/F"); subleadLepCostheta = 0;
    _outputTree->Branch("subleadLepRatio15", &subleadLepRatio15, "subleadLepRatio15/F"); subleadLepRatio15 = 0;
    _outputTree->Branch("subleadLepRatio30", &subleadLepRatio30, "subleadLepRatio30/F"); subleadLepRatio30 = 0;

    _outputTree->Branch("num_jet", &num_jet, "num_jet/I"); num_jet = 0;
    _outputTree->Branch("JetsInvMass", &JetsInvMass, "JetsInvMass/F"); JetsInvMass = 0;
    _outputTree->Branch("JetsRecoilMass", &JetsRecoilMass, "JetsRecoilMass/F"); JetsRecoilMass = 0;
    _outputTree->Branch("jet1cosTheta", &jet1cosTheta, "jet1cosTheta/F"); jet1cosTheta = 999;
    _outputTree->Branch("jet2cosTheta", &jet2cosTheta, "jet2cosTheta/F"); jet2cosTheta = 999;
    _outputTree->Branch("jet1En", &jet1En, "jet1En/F"); jet1En = 0;
    _outputTree->Branch("jet2En", &jet2En, "jet2En/F"); jet2En = 0;
    _outputTree->Branch("jet1MCPDG", &jet1MCPDG, "jet1MCPDG/I"); jet1MCPDG = 999;
    _outputTree->Branch("jet2MCPDG", &jet2MCPDG, "jet2MCPDG/I"); jet2MCPDG = 999;



    for (int i = 0; i < 5; ++i)
    {
        leadLepEn_501[i] = 0;
        leadLepEn_22[i] = 0;
        leadLepEn_21120[i] = 0;
        leadLepEn_other[i] = 0;

        subleadLepEn_501[i] = 0;
        subleadLepEn_22[i] = 0;
        subleadLepEn_21120[i] = 0;
        subleadLepEn_other[i] = 0;
    }
    _outputTree->Branch("leadLepEn_501", &leadLepEn_501, "leadLepEn_501[5]/F");
    _outputTree->Branch("leadLepEn_22", &leadLepEn_22, "leadLepEn_22[5]/F");
    _outputTree->Branch("leadLepEn_21120", &leadLepEn_21120, "leadLepEn_21120[5]/F");
    _outputTree->Branch("leadLepEn_other", &leadLepEn_other, "leadLepEn_other[5]/F");

    _outputTree->Branch("subleadLepEn_501", &subleadLepEn_501, "subleadLepEn_501[5]/F");
    _outputTree->Branch("subleadLepEn_22", &subleadLepEn_22, "subleadLepEn_22[5]/F");
    _outputTree->Branch("subleadLepEn_21120", &subleadLepEn_21120, "subleadLepEn_21120[5]/F");
    _outputTree->Branch("subleadLepEn_other", &subleadLepEn_other, "subleadLepEn_other[5]/F");


    leadConePDG = 999;
    leadConeSubEn = 0;
    leadConeAngle = 999;
    leadConeTotalMass = 0;
    leadConeTotalRatio = 0;
    _outputTree->Branch("leadConePDG", &leadConePDG, "leadConePDG/I");
    _outputTree->Branch("leadConeSubEn", &leadConeSubEn, "leadConeSubEn/F");
    _outputTree->Branch("leadConeAngle", &leadConeAngle, "leadConeAngle/F");
    _outputTree->Branch("leadConeTotalMass", &leadConeTotalMass, "leadConeTotalMass/F");
    _outputTree->Branch("leadConeTotalRatio", &leadConeTotalRatio, "leadConeTotalRatio/F");


    sametypePairMass = 0;
    difftypePairMass = 0;
    sametypePairRecoilMass = 0;
    difftypePairRecoilMass = 0;
    _outputTree->Branch("sametypePairMass", &sametypePairMass, "sametypePairMass/F"); 
    _outputTree->Branch("difftypePairMass", &difftypePairMass, "difftypePairMass/F");
    _outputTree->Branch("sametypePairRecoilMass", &sametypePairRecoilMass, "sametypePairRecoilMass/F"); 
    _outputTree->Branch("difftypePairRecoilMass", &difftypePairRecoilMass, "difftypePairRecoilMass/F");

    tauDecay = 999;
    _outputTree->Branch("tauDecay", &tauDecay, "tauDecay/I");

    for (int i = 0; i < 2; i++)
    {
        HbbL[i] = 999;
        HccL[i] = 999;
    }
    _outputTree->Branch("HbbL", &HbbL, "HbbL[2]/F");
    _outputTree->Branch("HccL", &HccL, "HccL[2]/F");


    Y12 = 999, Y23 = 999, Y34 = 999;
    _outputTree->Branch("Y12", &Y12, "Y12/F");
    _outputTree->Branch("Y23", &Y23, "Y23/F");
    _outputTree->Branch("Y34", &Y34, "Y34/F");

    for (int i = 0; i < 4; i++)
    {
        jet14m[i] = 999;
        jet24m[i] = 999;
    }
    _outputTree->Branch("jet14m", &jet14m, "jet14m[4]/F");
    _outputTree->Branch("jet24m", &jet24m, "jet24m[4]/F");


    // quark 1 and quark 2
    num_quark = 0;
    angle1 = 999;
    angle2 = 999;
    for (int i = 0; i < 4; i++)
    {
        quark14m[i] = 999;
        quark24m[i] = 999;
    }
    quark1PDG = 999;
    quark2PDG = 999;

    _outputTree->Branch("num_quark", &num_quark, "num_quark/I");
    _outputTree->Branch("angle1", &angle1, "angle1/F");
    _outputTree->Branch("angle2", &angle2, "angle2/F");
    _outputTree->Branch("quark14m", &quark14m, "quark14m[4]/F");
    _outputTree->Branch("quark24m", &quark24m, "quark24m[4]/F");
    _outputTree->Branch("quark1PDG", &quark1PDG, "quark1PDG/I");
    _outputTree->Branch("quark2PDG", &quark2PDG, "quark2PDG/I");

    // quark 1 and quark 2
    num_lepton = 0;
    for (int i = 0; i < 4; i++)
    {
        lepton14m[i] = 999;
        lepton24m[i] = 999;
    }
    lepton1PDG = 999;
    lepton2PDG = 999;
    _outputTree->Branch("num_lepton", &num_lepton, "num_lepton/I");
    _outputTree->Branch("lepton14m", &lepton14m, "lepton14m[4]/F");
    _outputTree->Branch("lepton24m", &lepton24m, "lepton24m[4]/F");
    _outputTree->Branch("lepton1PDG", &lepton1PDG, "lepton1PDG/I");
    _outputTree->Branch("lepton2PDG", &lepton2PDG, "lepton2PDG/I");


}
void vcb2::init()
{

    printParameters();
    TFile *tree_file = new TFile(_treeFileName.c_str(), (_overwrite ? "RECREATE" : "UPDATE"));
    if (!tree_file->IsOpen())
    {
        delete tree_file;
        tree_file = new TFile(_treeFileName.c_str(), "NEW");
    }

    _outputTree = new TTree(_treeName.c_str(), _treeName.c_str());    
    _outputTree->SetAutoSave(32 * 1024 * 1024); // autosave every 32MB

    setBranchAndValue(_outputTree);

    Num = 0;
}

TLorentzVector getV4(ReconstructedParticle *pfo)
{
    if (pfo)
    {
        return TLorentzVector(pfo->getMomentum(), pfo->getEnergy());
    }
    else
    {
        return TLorentzVector();
    }
}

#if 0
void updateSubLeadTrack(std::vector<ReconstructedParticle *> const &tracks,
 ReconstructedParticle * const leadTrack,
 ReconstructedParticle *&subTrack) {
    if(leadTrack) {
        for(size_t i = 0; i < tracks.size(); ++i) {

            ReconstructedParticle *candidate = tracks[i];
            if(leadTrack != candidate && candidate->getCharge() != leadTrack->getCharge()) {
                if(!subTrack) {
                    subTrack= candidate;
                } else if(subTrack->getEnergy() < candidate->getEnergy()){
                    subTrack = candidate;
                }
            }
        }
    }
}
//updateSubLeadTrack(vMuon, leadMuon, subTrack);
//updateSubLeadTrack(vElec, leadMuon, subTrack);

#endif

bool goodParticle(ReconstructedParticle *pfo) {
    //int pdgid = abs(pfo->getType());
    //if(pdgid == 11) return true;
    //if(pdgid == 13) return true;
    //if(pdgid == 22) return true;

    if(pfo->getEnergy() > 0.5) {
        return true;
    }

    //if(pdgid == 21120) return false;
    //if(pdgid == 501) return false;
    //if(pdgid == 21) return false;
    //if(pdgid == 31) return false;
    //if(pdgid == 20) return false;
    //if(pdgid == 101) return false;
    //if(pdgid == 1) return false;

    return false;
}

void vcb2::fillParticles(LCCollection *col_PFO,
    std::vector<ReconstructedParticle *> &vElec,
    std::vector<ReconstructedParticle *> &vMuon,
    std::vector<ReconstructedParticle *> &vPion,
    std::vector<ReconstructedParticle *> &vKaon,
    std::vector<ReconstructedParticle *> &vProton,
    std::vector<ReconstructedParticle *> &vGamma,
    TLorentzVector &TLPFO
    ) {

    int nPFO = col_PFO->getNumberOfElements();
    for (int i = 0; i < nPFO; i++)
    {
        ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(col_PFO->getElementAt(i));
        TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
        TLPFO += temp;
        int type = abs(pfo->getType());        

        //cout << "charge         " << pfo->getCharge() << endl;
        //cout << "pid    " << pfo->getType() << endl;
        //cout << "energy " << pfo->getEnergy() << endl;
        if(pfo->getCharge()) {
            nChargedPFOs += 1;
            if(goodParticle(pfo)) {
                nGoodChargedPFOs += 1;
            }
        }
        if(goodParticle(pfo)) {
            nGoodPFOs += 1;
        }

        if (type == 11)
        {
            vElec.push_back(pfo);
        }
        if (type == 13)
        {
            vMuon.push_back(pfo);
        }
        if (type == 211)
        {
            vPion.push_back(pfo);
        }
        if (type == 321)
        {
            vKaon.push_back(pfo);
        }
        if (type == 2212)
        {
            vProton.push_back(pfo);
        }
        if (type == 22)
        {
            vGamma.push_back(pfo);
        }
    }
    //if(nPFO > 10){
        //cout << "nChargedPFOs "  << nChargedPFOs << endl;
        //cout << "nGoodChargedPFOs " << nGoodChargedPFOs << endl;
        //cout << "nGoodPFOs " << nGoodPFOs << endl;
        //cout << "nPFO " << nPFO << endl;
        //getchar();
    //}

}


bool isQuark(int pdg) {
    return (abs(pdg) >= 1) && (abs(pdg) <= 8);
}
bool isLepton(int pdg) {
    return (abs(pdg) >= 11) && (abs(pdg) <= 18);
}

bool fromQuarkDecay(MCParticle *mcp) {

    for(;mcp;){
        if(isQuark(mcp->getPDG())) {
            return true;
        }
        if(mcp->getParents().size() == 0) {
            return false;
        }
        mcp = mcp->getParents()[0];
    }
    return false;
}

// how (first) primary tau decay
// 999: no tau find
// 1: hadronic decay
// 11,-11,13,-13: decay to electrons and muons
void TauDecayPDG(LCCollection *col_MCP, int & tauDecay)
{
    int n_MCP = col_MCP->getNumberOfElements();

    MCParticle *tauCandidate = NULL;
    for (int i = 0; i < n_MCP; i++)
    {
        MCParticle *a_MCP = dynamic_cast<MCParticle *>(col_MCP->getElementAt(i));
        int PDG = a_MCP->getPDG();

        if (abs(PDG) == 15 && !fromQuarkDecay(a_MCP)) // tauon
        {
            //cout << "tau find" << endl;
            for (int j = 0; j < (int)a_MCP->getDaughters().size(); j++)
            {
                MCParticle *tauDau = a_MCP->getDaughters()[j];
                if (abs(tauDau->getPDG()) == 16) // tau neutrino
                {
                    tauCandidate = a_MCP;
                    break;
                }
            }
        }
        if(tauCandidate) {
            break;
        }
    }

    if (tauCandidate != NULL)
    {
        //cout << "tauCandidate" << endl;
        bool hadr = true;
        for (int i = 0; i < (int)tauCandidate->getDaughters().size(); i++)
        {
            MCParticle *tauD = tauCandidate->getDaughters()[i];
            if (abs(tauD->getPDG()) == 11 || abs(tauD->getPDG()) == 13)
            {
                tauDecay = tauD->getPDG();
                hadr = false;
            }
        }
        if(hadr) {
            tauDecay = 1;
        }
    }

}

void vcb2::fillLeptons(std::vector<MCParticle*> &leptonvec) {
    // store quarks if there are two quarks
    // ordered by jets, if there are two jets
    num_lepton = leptonvec.size();
    cout << "num_lepton : " << num_lepton << endl;
    if (leptonvec.size() == 2)
    {
        MCParticle *quark1 = leptonvec.at(0);
        MCParticle *quark2 = leptonvec.at(1);

        TLorentzVector TLquark1(quark1->getMomentum(), quark1->getEnergy());
        TLorentzVector TLquark2(quark2->getMomentum(), quark2->getEnergy());
        TVector3 TVquark1 = TLquark1.Vect();
        TVector3 TVquark2 = TLquark2.Vect();

        lepton14m[0] = quark1->getMomentum()[0];
        lepton14m[1] = quark1->getMomentum()[1];
        lepton14m[2] = quark1->getMomentum()[2];
        lepton14m[3] = quark1->getEnergy();
        lepton24m[0] = quark2->getMomentum()[0];
        lepton24m[1] = quark2->getMomentum()[1];
        lepton24m[2] = quark2->getMomentum()[2];
        lepton24m[3] = quark2->getEnergy();

        lepton1PDG = quark1->getPDG();
        lepton2PDG = quark2->getPDG();

    }
}
void vcb2::fillJets(LCCollection *col_Jet,
std::vector<MCParticle*> &quarkvec)
{
    if(col_Jet) {
        num_jet = col_Jet->getNumberOfElements();

        TLorentzVector TLJets(0, 0, 0, 0);
        TLorentzVector TLCM(0, 0, 0, _centerOfMassEnergy);

        for (int i = 0; i < num_jet; ++i)
        {
            TLJets += lorentzV4(dynamic_cast<ReconstructedParticle *>(col_Jet->getElementAt(i)));
        }
        JetsInvMass = TLJets.M();
        JetsRecoilMass = (TLCM - TLJets).M();

        if (num_jet >= 1)
        {

            PIDHandler pidh(col_Jet);
            int algo = 0;

            algo = pidh.getAlgorithmID("yth");
            int iy12 = -1, iy23 = -1, iy34 = -1;
            iy12 = pidh.getParameterIndex(algo, "y12");
            iy23 = pidh.getParameterIndex(algo, "y23");
            iy34 = pidh.getParameterIndex(algo, "y34");
            ReconstructedParticle *jet = dynamic_cast<ReconstructedParticle *>(col_Jet->getElementAt(0));
            const ParticleID &pid = pidh.getParticleID(jet, algo);
            Y12 = pid.getParameters()[iy12];
            Y23 = pid.getParameters()[iy23];
            Y34 = pid.getParameters()[iy34];

            cout << "Ys: " << Y12 << " " << Y23 << " " << Y34 << endl;
        }

        if (num_jet == 2)
        {

            PIDHandler pidh(col_Jet);
            int algo = 0;
            int ibtag = -1;
            int ictag = -1;
            // int ibctag = -1;
            int icat = -1;
            algo = pidh.getAlgorithmID("lcfiplus");
            ibtag = pidh.getParameterIndex(algo, "BTag");
            ictag = pidh.getParameterIndex(algo, "CTag");
            // ibctag = pidh.getParameterIndex(algo, "BCTag");
            icat = pidh.getParameterIndex(algo, "Category");

            ReconstructedParticle *jet1 = dynamic_cast<ReconstructedParticle *>(col_Jet->getElementAt(0));

            const ParticleID &pid1 = pidh.getParticleID(jet1, algo);
            double btag1 = pid1.getParameters()[ibtag];
            double ctag1 = pid1.getParameters()[ictag];

            TLorentzVector TLjet1(jet1->getMomentum(), jet1->getEnergy());
            TVector3 TVjet1 = TLjet1.Vect();
            jet1cosTheta = TVjet1.CosTheta();
            jet1En = jet1->getEnergy();
            jet14m[0] = jet1->getMomentum()[0];
            jet14m[1] = jet1->getMomentum()[1];
            jet14m[2] = jet1->getMomentum()[2];
            jet14m[3] = jet1->getEnergy();

            ReconstructedParticle *jet2 = dynamic_cast<ReconstructedParticle *>(col_Jet->getElementAt(1));

            const ParticleID &pid2 = pidh.getParticleID(jet2, algo);
            double btag2 = pid2.getParameters()[ibtag];
            double ctag2 = pid2.getParameters()[ictag];

            TLorentzVector TLjet2(jet2->getMomentum(), jet2->getEnergy());
            TVector3 TVjet2 = TLjet2.Vect();
            jet2cosTheta = TVjet2.CosTheta();
            jet2En = jet2->getEnergy();

            jet24m[0] = jet2->getMomentum()[0];
            jet24m[1] = jet2->getMomentum()[1];
            jet24m[2] = jet2->getMomentum()[2];
            jet24m[3] = jet2->getEnergy();

            cout << "btag1 : " << btag1 << " btag2 " << btag2 << endl;
            cout << "ctag1 : " << ctag1 << " ctag2 " << ctag2 << endl;
            cout << "jet1cosTheta : " << jet1cosTheta << " " << jet1En << endl;
            cout << "jet2cosTheta : " << jet2cosTheta << " " << jet2En << endl;

            HbbL[0] = btag1;
            HccL[0] = ctag1;
            HbbL[1] = btag2;
            HccL[1] = ctag2;
        }
    }

    // store quarks if there are two quarks
    // ordered by jets, if there are two jets
    num_quark = quarkvec.size();
    cout << "num_quark : " << num_quark << endl;
    if (quarkvec.size() == 2)
    {
        MCParticle *quark1 = quarkvec.at(0);
        MCParticle *quark2 = quarkvec.at(1);

        TLorentzVector TLquark1(quark1->getMomentum(), quark1->getEnergy());
        TLorentzVector TLquark2(quark2->getMomentum(), quark2->getEnergy());
        TVector3 TVquark1 = TLquark1.Vect();
        TVector3 TVquark2 = TLquark2.Vect();


        if(num_jet == 2) {

            ReconstructedParticle *jet1 = dynamic_cast<ReconstructedParticle *>(col_Jet->getElementAt(0));
            ReconstructedParticle *jet2 = dynamic_cast<ReconstructedParticle *>(col_Jet->getElementAt(1));
            TLorentzVector TLjet1(jet1->getMomentum(), jet1->getEnergy());
            TVector3 TVjet1 = TLjet1.Vect();
            TLorentzVector TLjet2(jet2->getMomentum(), jet2->getEnergy());
            TVector3 TVjet2 = TLjet2.Vect();

            if (TVjet1.Angle(TVquark1) + TVjet2.Angle(TVquark2) > TVjet1.Angle(TVquark2) + TVjet2.Angle(TVquark1))
            {
                std::swap(quark1, quark2);
                std::swap(TLquark1, TLquark2);
                std::swap(TVquark1, TVquark2);
            }

            angle1 = TVjet1.Angle(TVquark1);
            angle2 = TVjet2.Angle(TVquark2);
        }

        quark14m[0] = quark1->getMomentum()[0];
        quark14m[1] = quark1->getMomentum()[1];
        quark14m[2] = quark1->getMomentum()[2];
        quark14m[3] = quark1->getEnergy();
        quark24m[0] = quark2->getMomentum()[0];
        quark24m[1] = quark2->getMomentum()[1];
        quark24m[2] = quark2->getMomentum()[2];
        quark24m[3] = quark2->getEnergy();

        quark1PDG = quark1->getPDG();
        quark2PDG = quark2->getPDG();

        jet1MCPDG = quark1PDG;
        jet2MCPDG = quark2PDG;
    }

}


void GetMCP(LCCollection *col_MCP, std::vector<MCParticle *> &quarkvec,
                                    std::vector<MCParticle *> &leptonvec)
{
    int n_MCP = col_MCP->getNumberOfElements();
    for (int i = 0; i < n_MCP; i++)
    {
        MCParticle *a_MCP = dynamic_cast<MCParticle *>(col_MCP->getElementAt(i));
        int NParents = a_MCP->getParents().size();
        int NDaughters = a_MCP->getDaughters().size();
        int PDG = a_MCP->getPDG();
        TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 tempV3 = temp.Vect();

        if (NParents == 0 && isQuark(PDG))
        {
            quarkvec.push_back(a_MCP);
        }
        if (NParents == 0 && isLepton(PDG))
        {
            leptonvec.push_back(a_MCP);
        }

        if (PDG == 25 && NDaughters == 4) // Higgs
        {
            // ...
        }
    }    
}


void vcb2::fracInCone(ReconstructedParticle *part, LCCollection *col_PFO, bool leading)
{
    int nPFO = col_PFO->getNumberOfElements();
    TLorentzVector seedMuon(part->getMomentum(), part->getEnergy());
    TVector3 seeMuonV3 = seedMuon.Vect();

    for (int i = 0; i < nPFO; i++)
    {
        ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(col_PFO->getElementAt(i));
        TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
        TVector3 TVtemp = temp.Vect();

        double angleInDegree = TVtemp.Angle(seeMuonV3)/3.1415926*180;        
        double energy = pfo->getEnergy();
        int type = pfo->getType();
        double angles[5] = {1,2,3,5,30};
        for(int i = 0; i < 5; ++i) {
            if(angleInDegree < angles[i]) {
                if(leading) {
                    if(type == 501) leadLepEn_501[i] += energy;
                    else if(type == 22) leadLepEn_22[i] += energy;
                    else if(type == 21120) leadLepEn_21120[i] += energy;
                    else if(type != _isoLepPDG) leadLepEn_other[i] += energy;
                } else {
                    if(type == 501) subleadLepEn_501[i] += energy;
                    else if(type == 22) subleadLepEn_22[i] += energy;
                    else if(type == 21120) subleadLepEn_21120[i] += energy;
                    else if(type != _isoLepPDG) subleadLepEn_other[i] += energy;
                }
            }
        }
    }
}

double energyInCone(ReconstructedParticle *part, LCCollection *col_PFO, double angleDeg)
{
    double enInCone = 0;
    int nPFO = col_PFO->getNumberOfElements();
    TLorentzVector seedMuon(part->getMomentum(), part->getEnergy());
    TVector3 seeMuonV3 = seedMuon.Vect();

    for (int i = 0; i < nPFO; i++)
    {
        ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(col_PFO->getElementAt(i));
        TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
        TVector3 TVtemp = temp.Vect();

        double angleInDegree = TVtemp.Angle(seeMuonV3)/3.1415926*180;        
        double energy = pfo->getEnergy();
        if (angleInDegree < angleDeg)
        {
            enInCone += energy;
        }
    }
    return enInCone;
}

double energyFraction(ReconstructedParticle *part, LCCollection *col_PFO, double angleDeg)
{
    return part->getEnergy() / energyInCone(part, col_PFO, angleDeg);
}

ReconstructedParticle *vcb2::saveIsoLepton(std::vector<ReconstructedParticle *> const &isoLeps, LCCollection *col_PFO)
{
    ReconstructedParticle *leadLep = NULL;
    if (isoLeps.size() > 0)
    {
        leadLep = isoLeps.at(0);

        leadLepM4[0] = leadLep->getMomentum()[0];
        leadLepM4[1] = leadLep->getMomentum()[1];
        leadLepM4[2] = leadLep->getMomentum()[2];
        leadLepM4[3] = leadLep->getEnergy();

        TLorentzVector leadLepV4 = TLorentzVector(leadLep->getMomentum(), leadLep->getEnergy());
        leadLepEn = leadLep->getEnergy();
        leadLepCharge = leadLep->getCharge();
        if(leadLep->getTracks().size() > 0) {
            Track *leadTrack = leadLep->getTracks()[0];
            leadD0 = leadTrack->getD0();
            leadZ0 = leadTrack->getZ0();
            leadIMP = pow(leadD0 * leadD0 + leadZ0 * leadZ0, 0.5);
        }
        TVector3 leadLepV3 = leadLepV4.Vect();
        leadLepCostheta = leadLepV3.CosTheta();

        leadLepRatio15 = energyFraction(leadLep, col_PFO, 15.);
        leadLepRatio30 = energyFraction(leadLep, col_PFO, 30.);
        fracInCone(leadLep, col_PFO, true);

        if (isoLeps.size() > 1)
        {
            ReconstructedParticle *subleadLep = isoLeps.at(1);

            subleadLepM4[0] = subleadLep->getMomentum()[0];
            subleadLepM4[1] = subleadLep->getMomentum()[1];
            subleadLepM4[2] = subleadLep->getMomentum()[2];
            subleadLepM4[3] = subleadLep->getEnergy();


            subleadLepEn = subleadLep->getEnergy();
            subleadLepCharge = subleadLep->getCharge();
            if(subleadLep->getTracks().size() > 0) {
                Track *subleadTrack = subleadLep->getTracks()[0];
                subleadD0 = subleadTrack->getD0();
                subleadZ0 = subleadTrack->getZ0();
                subleadIMP = pow(subleadD0 * subleadD0 + subleadZ0 * subleadZ0, 0.5);
            }
            subleadLepRatio15 = energyFraction(subleadLep, col_PFO, 15.);
            subleadLepRatio30 = energyFraction(subleadLep, col_PFO, 30.);
            fracInCone(subleadLep, col_PFO, false);
        }
    }
    return leadLep;
}


void vcb2::doProcessEvent(LCEvent *evtP)
{

    setBranchAndValue(NULL); // default values

    cout << "Next Event *******************************************************************************************************" << endl;
    eventNr = evtP->getEventNumber();
    cout << "eventNr : " << eventNr << " Num : " << Num << endl;

    TLorentzVector TLCM(0, 0, 0, _centerOfMassEnergy);

    LCCollection *col_PFO = evtP->getCollection("ArborPFOs");

    int nPFO = col_PFO->getNumberOfElements();

    TLorentzVector TLPFO(0, 0, 0, 0);
    std::vector<ReconstructedParticle *> vElec;
    std::vector<ReconstructedParticle *> vMuon;
    std::vector<ReconstructedParticle *> vPion;
    std::vector<ReconstructedParticle *> vKaon;
    std::vector<ReconstructedParticle *> vProton;
    std::vector<ReconstructedParticle *> vGamma;

    fillParticles(col_PFO, vElec, vMuon, vPion, vKaon, vProton, vGamma, TLPFO);

    TLorentzVector missV4 = TLCM - TLPFO;
    visEn = TLPFO.E();
    missPt = (missV4.Vect()).Perp();
    missM = missV4.M();

    multiplicity = nPFO;
    nPFOs = nPFO;

    cout << "nChargedPFOs "  << nChargedPFOs << endl;
    cout << "nGoodChargedPFOs " << nGoodChargedPFOs << endl;
    cout << "nGoodPFOs " << nGoodPFOs << endl;
    cout << "nPFO " << nPFO << endl;

    sort(vElec.begin(), vElec.end(), sortEn);
    sort(vMuon.begin(), vMuon.end(), sortEn);
    sort(vPion.begin(), vPion.end(), sortEn);
    sort(vKaon.begin(), vKaon.end(), sortEn);
    sort(vProton.begin(), vProton.end(), sortEn);
    sort(vGamma.begin(), vGamma.end(), sortEn);

    cout << "the total energy of PFOs is : " << TLPFO.E() << endl;
    cout << "vElec.size() : " << vElec.size() << endl;
    cout << "vMuon.size() : " << vMuon.size() << endl;
    cout << "vPion.size() : " << vPion.size() << endl;
    cout << "vKaon.size() : " << vKaon.size() << endl;
    cout << "vProton.size() : " << vProton.size() << endl;
    cout << "vGamma.size() : " << vGamma.size() << endl;

    if (vMuon.size() > 0)
    {
        ReconstructedParticle *leadMuon = vMuon.at(0);
        leadMuonEn = leadMuon->getEnergy();
    }
    if (vElec.size() > 0)
    {
        ReconstructedParticle *leadElec = vElec.at(0);
        leadElecEn = leadElec->getEnergy();
    }
    if (vPion.size() > 0)
    {
        ReconstructedParticle *leadPion = vPion.at(0);
        leadPionEn = leadPion->getEnergy();
    }
    if (vGamma.size() > 0)
    {
        ReconstructedParticle *leadPion = vGamma.at(0);
        leadGammaEn = leadPion->getEnergy();
    }

    cout << "leadElecEn " << leadElecEn << endl;
    cout << "leadMuonEn " << leadMuonEn << endl;

    ReconstructedParticle *leadLep = NULL;
    if (abs(_isoLepPDG) == 13)
        leadLep = saveIsoLepton(vMuon, col_PFO);
    else if (abs(_isoLepPDG) == 11)
        leadLep = saveIsoLepton(vElec, col_PFO);

    std::vector<ReconstructedParticle*> *homoset = &vElec;
    std::vector<ReconstructedParticle*> *heteset = &vMuon;
    std::vector<ReconstructedParticle*> vElecMuon;
    vElecMuon.insert(vElecMuon.end(), vElec.begin(), vElec.end());
    vElecMuon.insert(vElecMuon.end(), vMuon.begin(), vMuon.end());

    if(abs(_isoLepPDG) == 13) {
        std::swap(homoset, heteset);
    }
        
    if(leadLep) {
        ReconstructedParticle* subleadCone = NULL;
        double maxEnergy = 0;
        TVector3 center = lorentzV4(leadLep).Vect();
        //for(int i = 0; i < (int)vElecMuon.size(); ++i) {
        for(int i = 0; i < (int)col_PFO->getNumberOfElements(); ++i) {

            //ReconstructedParticle* part = vElecMuon.at(i);
            ReconstructedParticle* part = 
                dynamic_cast<ReconstructedParticle*>( col_PFO->getElementAt(i) );

            if(part->getCharge() != 0) {
                if(center.Angle(lorentzV4(part).Vect()) < 30.0/180.0*3.1415926) {
                    if(part == leadLep) continue;
                    if(part->getEnergy() > maxEnergy) {
                        maxEnergy = part->getEnergy();
                        subleadCone =  part;
                    }
                }
            }
        }

        //if(subleadCone && (subleadCone->getCharge() != leadLep->getCharge())
        //    && (abs(subleadCone->getType()) == abs(leadLep->getType())) ) {
        cout << "subleadCone " << subleadCone << endl;
        if(subleadCone && (subleadCone->getCharge() != leadLep->getCharge())) {

            double total_energy = energyInCone(leadLep, col_PFO, 30.);
            cout << "total energy " <<  total_energy << endl;
            leadConeSubEn = subleadCone->getEnergy();
            leadConeAngle = lorentzV4(subleadCone).Vect().Angle(center);
            leadConePDG = subleadCone->getType();
            leadConeTotalMass = (lorentzV4(subleadCone) + lorentzV4(leadLep)).M();
            leadConeTotalRatio = (subleadCone->getEnergy() + leadLep->getEnergy()) / total_energy;
        }



        if((int)homoset->size() >= 2) {
            ReconstructedParticle* sublead = homoset->at(1);
            if(sublead->getCharge() != leadLep->getCharge()) {
                sametypePairMass = (lorentzV4(sublead) + lorentzV4(leadLep)).M();
                sametypePairRecoilMass = (TLCM - lorentzV4(sublead) - lorentzV4(leadLep)).M();
            }
        }
        if((int)heteset->size() >= 1) {
            ReconstructedParticle* sublead = heteset->at(0);
            if(sublead->getCharge() != leadLep->getCharge()) {
                difftypePairMass = (lorentzV4(sublead) + lorentzV4(leadLep)).M();
                difftypePairRecoilMass = (TLCM - lorentzV4(sublead) - lorentzV4(leadLep)).M();
            }
        }
    }

    
    LCCollection *col_MCP = evtP->getCollection("MCParticle");
    std::vector<MCParticle *> quarkvec;
    std::vector<MCParticle *> leptonvec;
    GetMCP(col_MCP, quarkvec, leptonvec);

    TauDecayPDG(col_MCP, tauDecay);
    cout << "tauDecay: " << tauDecay << endl;

    LCCollection *col_Jet = evtP->getCollection("RefinedJets");
    fillJets(col_Jet, quarkvec);
    fillLeptons(leptonvec);
}

void vcb2::processEvent(LCEvent *evtP)
{

    if (evtP)
    {
        try {
            doProcessEvent(evtP);
        }
        catch (lcio::DataNotAvailableException err)
        {
        }
    }
    // fill anyway!
    _outputTree->Fill();
    Num++;
    cout << " ***** process event end ****" << endl;
}

void vcb2::end()
{

    if (_outputTree)
    {

        TFile *tree_file = _outputTree->GetCurrentFile(); // just in case we switched to a new file
        // tree_file->cd();
        tree_file->Write();
        delete tree_file;
        // tree_file->Close();
    }
}
