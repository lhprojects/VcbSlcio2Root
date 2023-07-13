
#include "StdHepReader.h"
#include "marlin/ProcessorMgr.h"

#include "IMPL/LCEventImpl.h"
#include "IMPL/LCRunHeaderImpl.h"

#include "UTIL/LCStdHepRdr.h"
#include "UTIL/LCTOOLS.h"


#include <iostream>
#include <vector>
#include <cctype>

std::string trim(const std::string& str) {
    size_t firstNonSpace = str.find_first_not_of(" \t\n\r");
    size_t lastNonSpace = str.find_last_not_of(" \t\n\r");
    
    if (firstNonSpace == std::string::npos || lastNonSpace == std::string::npos) {
        // The string is empty or contains only whitespace
        return "";
    }
    
    return str.substr(firstNonSpace, lastNonSpace - firstNonSpace + 1);
}



StdHepReader2 aStdHepReader ;
StdHepReader2::StdHepReader2() : DataSourceProcessor("StdHepReader2") {
  
  _description = "Reads StdHep files as input and creates LCIO events with MCParticle collections."
    " Make sure to not specify any LCIOInputFiles in the steering in order to read StdHep files." ;
  
  registerProcessorParameter("StdHepFileName" , 
      "input files"  ,
      _fileNames,
      _fileNames) ;
  
}

StdHepReader2*  StdHepReader2::newProcessor() { 
  return new StdHepReader2 ;
}

void StdHepReader2::init() {    
  printParameters() ;    
}


void StdHepReader2::readDataSource( int numEvents ) {
  
  for(int i = 0; i < (int)_fileNames.size(); ++i) {

      std::string filename = _fileNames[i];

        LCStdHepRdr *rdr = new LCStdHepRdr(filename.c_str());

        LCCollection *col;
        LCEventImpl *evt;

        int evtNum = 0;
        int runNum = 0;

        while ((col = rdr->readEvent()) != 0)
        {

          if (numEvents > 0 && evtNum + 1 > numEvents)
          {
            delete col;
            break;
          }

          if (isFirstEvent())
          { // create run header

            LCRunHeaderImpl *rHdr = new LCRunHeaderImpl;

            rHdr->setDescription(" Events read from stdhep input file: " + filename);
            rHdr->setRunNumber(runNum);

            ProcessorMgr::instance()->processRunHeader(rHdr);
            _isFirstEvent = false;
          }

          evt = new LCEventImpl;
          evt->setRunNumber(runNum);
          evt->setEventNumber(evtNum++);

          evt->addCollection(col, "MCParticle");

          ProcessorMgr::instance()->processEvent(evt);

          delete evt;
        }

        delete rdr;
    }

}

void StdHepReader2::end() {

}


