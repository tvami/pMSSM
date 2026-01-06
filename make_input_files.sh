#!/bin/bash

CHUNK_SIZE=10000

# Define datasets and corresponding output files
declare -A DATASETS=(
    ["/PMSSM_set_2_LL_1_TuneCP2_13TeV-pythia8/RunIISummer20UL17NanoAODv9-FSMiniUL17_NANOv9_FSUL17_106X_mc2017_realistic_v9_ext1-v2/NANOAODSIM"]="inputdataset_PMSSM_set_2_LL_1_2017.txt"
    ["/PMSSM_set_2_LL_2_TuneCP2_13TeV-pythia8/RunIISummer20UL17NanoAODv9-FSMiniUL17_NANOv9_FSUL17_106X_mc2017_realistic_v9_ext1-v2/NANOAODSIM"]="inputdataset_PMSSM_set_2_LL_2_2017.txt"
    ["/PMSSM_set_1_LL_TuneCP2_13TeV-pythia8/RunIISummer20UL18NanoAODv9-FSMiniUL18_NANOv9_FSUL18_106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM"]="inputdataset_PMSSM_set_1_LL_2018.txt"
    ["/PMSSM_set_2_LL_1_TuneCP2_13TeV-pythia8/RunIISummer20UL18NanoAODv9-FSMiniUL18_NANOv9_FSUL18_106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM"]="inputdataset_PMSSM_set_2_LL_1_2018.txt"
    ["/PMSSM_set_2_LL_2_TuneCP2_13TeV-pythia8/RunIISummer20UL18NanoAODv9-FSMiniUL18_NANOv9_FSUL18_106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM"]="inputdataset_PMSSM_set_2_LL_2_2018.txt"
)

# Process each dataset
for DATASET in "${!DATASETS[@]}"; do
    OUTPUT="${DATASETS[$DATASET]}"
    
    echo "Querying files for dataset: $DATASET"
    echo "Writing to: $OUTPUT"
    echo ""
    
    # Clear output file if it exists
    > $OUTPUT
    
    # Get list of files
    FILES=$(dasgoclient -query="file dataset=$DATASET")
    
    # Loop through each file and get the number of events
    for FILE in $FILES; do
        NEVENTS=$(dasgoclient -query="file=$FILE" -json | grep -o '"nevents":[0-9]*' | cut -d':' -f2)
        echo "Processing: $FILE ($NEVENTS events)"
        
        # Break into 10k chunks
        START=0
        while [ $START -lt $NEVENTS ]; do
            END=$((START + CHUNK_SIZE))
            if [ $END -gt $NEVENTS ]; then
                END=$NEVENTS
            fi
            echo "$START $END $FILE" >> $OUTPUT
            START=$END
        done
    done
    
    echo "Done with $OUTPUT"
    echo ""
done

echo "All datasets processed!"