# Bipolar Diagnostic Test Replication
A reliable diagnostic genetic test can potentially aid early and accurate detection of bipolar disorder and reduce delays to treatment, which are currently high. 

A previous study showed promising results for using four mutations (DiO1, DiO2 Gly3Asp, DiO2 Thr92Ala and SLCO1C1) as diagnostic biomarkers for bipolar disorder. This study compared the prevalence of these four mutations in n=199 patients with bipolar disorder with an open-source healthy control sample and a cohort of n=179 patients with depression (for DiO2 Thr92Ala only). Three of these mutations (DiO1, DiO2 Thr92Ala and SLCO1C1) were significantly more prevalent in patients with bipolar disorder compared to healthy controls and the prevalence of DiO2 Thr92Ala was similarly higher in patients with bipolar disorder compared to those with depression. These mutations were able to discriminate between patients with bipolar disorder and healthy controls with high sensitivity (up to 86.7%) and acceptable specificity (up to 65.7%), with similar results discriminating between patients with bipolar disorder and depression. The overall balanced accuracy was not reported.

Given the paucity of replication research in psychiatry, this study primarily aims to test the generalisability of these findings in a larger sample from the UK Biobank Resource. This step is an essential requirement to validate potential diagnostic genetic tests for bipolar disorders. 

This script:
1. Reads in raw phenotypic and genetic data from the UK Biobank Resource
2. Performs data processing and cleaning
3. Generates a sociodemographics table
4. Compares prevalence of mutations across groups
5. Calculates balanced accuracy, sensitivity, specificity, positive predictive values and negative predictive values for each individual mutation, any mutation and all mutations in discriminating between patients with bipolar disorder and i) patients with depression; ii) healthy controls
