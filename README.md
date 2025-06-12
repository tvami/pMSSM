
 # Make the input eff templates with the LLP analzyer using the ntuples
```
make
./RazorRun lists/short.txt  HSCPAnalyzer_pMSSM -f=PPStau_M-557_new3.root -d=no
./RazorRun lists/localStauList.txt HSCPAnalyzer_pMSSM -f=triggerAndPres.root  -d=no
```

 # Run on NANOAOD 
```
root -l -q 'analyze_pmssm.C("PMSSM_set_1_LL_TuneCP2_13TeV_FS_NANOv9.root")'
```
