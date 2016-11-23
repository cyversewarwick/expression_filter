# Expression Filter

## The Purpose of the App

On the whole, CyVerse tools working with time series expression data accept data in a similarly formatted fashion. However, you may want to supply different parts of the data to different apps - some require both control and treated, with as many replicates as possible, while others want just differentially expressed treated genes averaged to a single value per time point. The documentation for each individual app features the details on what expression data it expects on input.

This helper app has been created to allow for easy reformatting of the data, depending on each app's desires. The operations performed here could be easily done manually in Excel or some other platform, but the app is nevertheless provided for convenience (and ease of chaining operations in workflows).

## Basic Input/Output

Expression Filter accepts an expression CSV on input and returns an expression CSV on output, with the contents slightly rejiggered according to the provided instructions.

## Test Run

In case you want to try out the different features, a demonstration expression CSV has been provided at `cyverseuk/expression_filter_testdata/input.csv` under Community Data. The folder also features an arbitrary gene list (`genes.txt`) and manipulated GP2S scores so that some genes come up as insignificant (`scores.txt`) if you want to try the relevant features.

## Input in Detail

### Expression CSV File

**Mandatory input.** In fact, the very reason why this app exists. First column for gene IDs, the first row may be condition information (with the relevant condition name repeating for each applicable column). If the condition information is provided, then the second row is to be information on time points. If condition information is absent, then the time point information is to be the first row.

### GP2S Scores

Output of the GP2S app ran with the same expression file. If provided, the expression data will be filtered to differentially expressed genes.

### Score Threshold

**Default:** 5

The GP2S score threshold to be used for differential expression identification. Genes with scores above this threshold are deemed differentially expressed. If no GP2S scores file is provided then this field is ignored by the app.

### Filter Gene List

One line per gene ID matching the ones provided in the expression CSV file. The file may be tab delimited and contain other information, but in that case the first column needs to be the relevant gene IDs with all irrelevant IDs removed. If supplied, the expression CSV file will be filtered down to the gene IDs present in this file. An example use would be filtering the expression data to just transcription factors.

### Remove Condition

Only provide this information if the first line of your expression file is treatment information. All data from the provided condition will be removed, and if the resulting expression file has only one condition remaining then the treatment line will be deleted. An example use would be removing control data from an expression file used to identify differentially expressed genes with GP2S in preparation for clustering with TCAP or BHC.

### Average Replicates

If checked, all the replicates for a given time point in a given condition will be averaged to a single value. Some tools (such as TCAP and Wigwams) expect their data to be formatted this way on input.