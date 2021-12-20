##run this command under each folder and the 10 trajectory RMSDs will be combined for easier analysis.

paste -d "," PMEMD{1..10}/RMSDresults.csv > combinedRMSD.csv
