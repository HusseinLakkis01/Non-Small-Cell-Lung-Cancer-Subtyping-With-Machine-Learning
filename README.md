# Non-Small-Cell-Lung-Cancer-Subtyping-With-Machine-Learning
Early diagnosis and prognosis of the type of cancer have become a necessity in cancer research, which can facilitate the subsequent clinical management of patients. Scientists have used different methods, such as early-stage screening, to detect types of cancer before they cause symptoms. Besides, they have developed new strategies for the early prediction of cancer treatment outcomes. With the advent of new technologies in the field of medicine, large amounts of data on cancer have been collected and are available to the medical research community.

Machine learning approaches are sufficient to assess a clinical diagnosis. Their final aim is to obtain trained models that predict the type and developmental character of malignant growth by using one or more classification attributes. Such algorithms can be used by doctors as an additional method to handle large amounts of patient data for diagnostic purposes. 

That being said, this code attempts to build accurate machine learning models that would classify NSCLCs into the two major subtypes.

## Datasets:
We used two publicly available data sets related to gene expression to test the efficacy of different machine learning algorithms. This data is obtained through Lung Cancer Explorer.

Lung Cancer Explorer (LCE) is an online tool developed by the Quantitative Biomedical Research Center (QBRC) of the UT Southwestern Medical Center that allows you to explore information on gene expression from hundreds of public lung cancer databases online for free

### 1. TCGA LUAD 2016, Cancer Genome Atlas Research Group: 
Comprehensive molecular profiling of lung adenocarcinoma published in Nature (https://www.ncbi.nlm.nih.gov/pubmed/?term=25079552). This dataset contains information from 576 patient tissue samples (517 tumors, and 59 normal). The data contains clinical information about the patients, in addition to the gene expression of each patient and some of their normal counterparts (20429 probes). The patients were segmented into their respective stage and substage and many factors were included.
### 2. TCGA LUSC 2016, Cancer Genome Atlas Research Group: 
Comprehensive molecular profiling of squamous cell lung cancer and published in Nature (https://www.ncbi.nlm.nih.gov/pubmed/?term=22960745). Similar to the previous set, and under the same conditions, this dataset contains both clinical and gene expression information from 552 patient tissue samples (501 tumors, and 51 normal).

## Procedure

### Minimal Redundancy, Maximal Relevance (MRMR):
Minimum Redundancy Feature Selection is an algorithm commonly used in a method to accurately identify the relationship between features and to narrow down their relevance to a specific outcome variable and is usually described as a Minimum Redundancy Maximum Relevance (MRMR).

Features can be chosen so that they are mutually distant from each other (not redundant) while still having a "high" correlation to the classification variable. This scheme, referred to as the Minimum Redundancy Maximum Relevance (MRMR) selection, was found to be more powerful than the maximum relevance selection.
The features are then ranked according to their score, i.e importance relative to the outcome, and top n features are selected. In this project, we used the package
‘praznik’, https://gitlab.com/mbq/praznik, which offers a function that allows parallelized execution for efficient performance.

### Principle Component Analysis:
PCA is very effective in reducing the dimensionality of data and thus can contribute to better, simpler, and more stable models. In this case, it would aid in the noise reduction of gene expression arrays. Data was then split into train and test sets.

### Model Training:
SVMs, Random Forests, KNN, ANN, and Naive Bayes models were built and trained on the principle component variables. Training and Test Errors were reported. Alsp. f1, precision, recall, and others were reported on the test set.
