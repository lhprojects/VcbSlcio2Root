/**
 * @brief Marlin processor for finding isolated leptons.
 * @author Ryo Yonamine <yonamine@post.kek.jp>
 * @author Tomohiko Tanabe <tomohiko@icepp.s.u-tokyo.ac.jp>
 * @version $Id:$
 *
 * Given a list of ReconstructedParticle, identify isolated leptons
 * based on the track cone energy, lepton identification,
 * and the track impact parameters (optional).
 */
#ifndef MCParticleFSParticle_h
#define MCParticleFSParticle_h 1

#include <string>
#include <vector>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

class MCParticleFSParticle : public Processor {

	public:

		virtual Processor*  newProcessor() { return new MCParticleFSParticle ; }

		MCParticleFSParticle() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;

	protected:

		/** Input collection */
		std::string _inputPFOsCollection;
        std::string _outmcpsimufsp;
        std::string _outMCPSIMURelation;
          int _replaceFL;
//        std::string _out941;
//        std::string _out942;
		LCCollection* _pfoCol;

} ;

#endif

