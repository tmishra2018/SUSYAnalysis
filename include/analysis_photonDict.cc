// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME analysis_photonDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "analysis_rawData.h"
#include "analysis_photon.h"
#include "analysis_ele.h"
#include "analysis_muon.h"
#include "analysis_scalefactor.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace {
  void TriggerDictionaryInitialization_analysis_photonDict_Impl() {
    static const char* headers[] = {
"analysis_rawData.h",
"analysis_photon.h",
"analysis_ele.h",
"analysis_muon.h",
"analysis_scalefactor.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/lcg/root/6.30.03-ca7ca986842b225f6fc22ae84d705ed8/include/",
"/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/include/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "analysis_photonDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "analysis_photonDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "analysis_rawData.h"
#include "analysis_photon.h"
#include "analysis_ele.h"
#include "analysis_muon.h"
#include "analysis_scalefactor.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("analysis_photonDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_analysis_photonDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_analysis_photonDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_analysis_photonDict() {
  TriggerDictionaryInitialization_analysis_photonDict_Impl();
}
