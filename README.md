# BMI3_isosplice
FInd and score the SJ using .bam file


you can use these to run the isosplice:
```
python main.py findSJ input.bam
python main.py buildDAG SJ_scores.csv
python main.py simplifyDAG DAG_edge.csv
python main.py scoreSJ cluster_scores.csv cluster_label.csv
```
