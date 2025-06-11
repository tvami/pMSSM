#ifndef DEF_HSCPAnalyzer_pMSSM
#define DEF_HSCPAnalyzer_pMSSM

#include "RazorAnalyzer.h"

class HSCPAnalyzer_pMSSM: public RazorAnalyzer {
    public: 
        HSCPAnalyzer_pMSSM(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
