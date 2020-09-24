# TUTORIAL FOR RGCCA R-SHINY

##### Version: 1.0

##### Author: Etienne CAMENEN

##### Key-words: 
omics, RGCCA, multi-block

##### EDAM operation: 
analysis, correlation, visualisation

##### Contact: 
arthur.tenenhaus@l2s.centralesupelec.fr

##### Short description
Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.) and produces textual and graphical outputs (e.g. variables and individuals plots).

---
## Contents
  - [Description](#description)
  - [1. Load the inputs ('Data' parameter tab)](#1-load-the-inputs-data-parameter-tab)
  - [2. Analysis parameters ('RGCCA' parameter tab)](#2-analysis-parameters-rgcca-parameter-tab)
    - [2.1. Analysis methods](#21-analysis-methods)
    - [2.2. Number of components and scaling](#22-number-of-components-and-scaling)
    - [2.3. Connection between blocks](#23-connection-between-blocks)
      - [2.3.1. Loading a connection file](#231-loading-a-connection-file)
      - [2.3.2. Superblock](#232-superblock)
      - [2.3.3. Supervised analysis](#233-supervised-analysis)
    - [2.4. Other R/SGCCA parameters](#24-other-rsgcca-parameters)
      - [2.4.1. Shrinkage parameter (Tau)](#241-shrinkage-parameter-tau)
      - [2.4.2. Sparsity coefficient](#242-sparsity-coefficient)
      - [2.4.3. Scheme function (advanced users)](#243-scheme-function-advanced-users)
    - [2.5. Pre & post-analysis functionalities](#25-pre--post-analysis-functionalities)
      - [2.5.1. Number of bootstraps, permutations, folds](#251-number-of-bootstraps-permutations-folds)
  - [3. Graphical parameters ('Graphic' parameter tab)](#3-graphical-parameters-graphic-parameter-tab)
    - [3.1. Display names](#31-display-names)
    - [3.2. Block (for the x/y-axis)](#32-block-for-the-xy-axis)
    - [3.3. Components (for the x/y-axis)](#33-components-for-the-xy-axis)
    - [3.4. Color the samples](#34-color-the-samples)
    - [3.5. Number of top variables](#35-number-of-top-variables)
    - [3.6. Save the graphics](#36-save-the-graphics)
    - [3.7. Dynamic graph parameters header](#37-dynamic-graph-parameters-header)
    - [3.8. Bootstrap indexes for x- and y- axis](#38-bootstrap-indexes-for-x--and-y--axis)
    - [3.9 Display cross-validation](#39-display-cross-validation)
  - [4. Visualize the plot tabs](#4-visualize-the-plot-tabs)
    - [4.1. Connection between blocks](#41-connection-between-blocks)
    - [4.2. Average variance explained (AVE)](#42-average-variance-explained-ave)
    - [4.3. Samples](#43-samples)
    - [4.4. Corcircle](#44-corcircle)
    - [4.5. Top variables](#45-top-variables)
    - [4.5. Bootstrap](#45-bootstrap)
    - [4.6. Permutation](#46-permutation)

## Description


We consider J data matrices X1 ,..., XJ. Each n × pj data matrix Xj = [ xj1, ..., xjpj ] is called a block and represents a set of pj variables observed on n individuals. The number and the nature of the variables may differ from one block to another, but the individuals must be the same across blocks. We assume that all variables are centered. The objective of RGCCA is to find, for each block, a weighted composite of variables (called block component) yj = Xj . aj, j = 1 ,..., J (where aj is a column-vector with pj elements) summarizing the relevant information between and within the blocks. The block components are obtained such that (i) block components explain well their own block and/or (ii) block components that are assumed to be connected are highly correlated. In addition, RGCCA integrates a variable selection procedure, called SGCCA, allowing the identification of the most relevant features (see [here](https://github.com/rgcca-factory/RGCCA#algorithm) for more information).


## 1. Load the inputs ('Data' parameter tab)

Before reading this tutorial, the Shiny application have to be [installed](https://github.com/rgcca-factory/RGCCA#installation) and [executed](https://github.com/rgcca-factory/RGCCA#execution).
The ```inst/extdata``` folder includes three blocks with the same individuals (corresponding to the countries here) but different types of variables (agriculture, industry and politic variables). In this dataset, according to Russett (1964), a high agriculture inequality and a low industrial development lead to an unstable political regime. Load the three working examples ```agriculture.tsv```, ```industry.tsv``` and ```politic.tsv``` in the ```blocks``` box (**Fig. 1**) (CTRL + click must be used for multiple file selection). The accepted format is one (for PCA) or multiple CSV files containing a matrix with:
- quantitative values only, with decimals separated by '.' and missing values labelled as "NA"
- samples in rows, labelled in the 1rst column with the same sample names between blocks (some samples could be missing in some blocks)
- variables in columns, labelled in the 1rst line without duplications in variable names between blocks

This format recommendation could be viewed with **the mouse over the question mark symbol** on the right of the file loading box.

By default, the character used in ```column separator``` parameter is the ```tabulation```. Change the separator to another one (e.g.,```semicolon```) to observe an error notification: 
> "politic block file has only one column. Check the separator."

![load_file_text](../../img/shiny/loadData.png)
 
*Fig. 1 : File loading panel (on the top right). The browsing box is used to load the blocks and the last one, to select the column separator.*


## 2. Analysis parameters ('RGCCA' parameter tab)

The analyse parameters are all set by default and the user could directly click on the ```run analysis``` button. To directly visualize the outputs, see the [last section](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#4-visualize-the-plot-tabs).

### 2.1. Analysis methods

After loading the data, a ```RGCCA``` tab will appear (**Fig. 2**). By default, the selected ```analysis method``` is set on ```RGCCA```. This tutorial will be focused on the RGCCA case, but another methods could be selected. When only one block file is loaded in the previous step, a ```PCA``` will be performed. By using two blocks, the interface will allow to select two-blocks analysis method (```PLS``` is selected by default). 

![analysis](../../img/shiny/analysis.png)

*Fig. 2 : The second parameter panel shows various options to customize the analysis: choose the analysis and the number of components, scale the blocks, choose a shrinkage, use the superblock or a supervised approach, choose a link function. In this example, the mouse over the question mark of scheme function parameter makes a help message appears*

### 2.2. Number of components and scaling

With all analysis methods, the ```number of components``` could be changed. By default, each blocks is set to two components (for 2D-graphics) or to an only component for an univariate block. The maximum of components allowed is limited by the minimum number of columns between all datasets. In the case of Russet data, two components are allowed because of the two columns in the industry block. In any cas, five components are the maximum value allowed. 
 
One could also selects ```scale [/ unscale] the blocks```. Either the option is selected or not, a data centering step is always performed. If selected, each block is normalized and then divided by the square root of its number of variables. When the data are already scaled, this step could be avoided by disabling the parameter.

### 2.3. Connection between blocks

This parameters are only accessible with R/SGCCA.

#### 2.3.1. Loading a connection file

The downloaded folder contains a design matrix (```connection.tsv```) corresponding to the relationship between each block: 1 if two blocks are connected and 0 otherwise. The expected format should be tabulation-separated and have column and row names corresponding to the blocks names. It is a symmetric matrix with the same dimension as the number of blocks. This file allows to add *a priori* information of correlation hypothesis between the blocks. It will not be taken in account with a superblock (see next section). After disabling the ```use a superblock``` option, load this file into the ```Connection design``` file box and observe the result on the plots. The ```connection.tsv``` file contains 1 in all non-diagonal cells and makes the assumption that all the blocks are related. 

A ```connection.tsv``` file will be automatically generated in the ```inst/shiny``` folder after running the analyse. This file corresponding to the connection matrix used by the model can be modified manually and reloaded in the software.

#### 2.3.2. Superblock 
By default, all the blocks are connected to a superblock, a concatenation of all the other blocks. The space spanned by global components is viewed as a compromise space that integrated all the modalities and facilitates the visualization of the results and their interpretation. To visualize the blocks without the superblock, disable the ```Use a superblock``` option.

#### 2.3.3. Supervised analysis
By selecting ```supervised analysis``` option, a drop-down menu appears to select the block used as a response. By selecting this block, all other blocks (predictors) will be only connected to this block. For example, select the ```agriculture``` block.

If a superblock is used, supervised analysis is automatically disabled, and inversely.

### 2.4. Other R/SGCCA parameters

#### 2.4.1. Shrinkage parameter (Tau)
By selecting a RGCCA,```use an optimal tau``` is disabled and each tau automatically set to 1 for each block (**Fig. 3**). One could make ```tau``` varying for each block from 1 (maximize the correlation between the variables of the selected block) to 0 (maximize the covariance). An optimal tau could not be set with a superblock configuration.
 
#### 2.4.2. Sparsity coefficient
By selecting a SGCCA, the ``` sparsity``` could be applied to each block. This coefficient varies from the inverse of the square root of the number of columns (the smaller set of variables) to 1 (all the variables are included).

Move again the cursor to an upper sparsity value to make it works.
 
#### 2.4.3. Scheme function (advanced users)
```Scheme function``` allows to select the link (i.e. scheme) function for covariance maximizations between block components among: 
- identity (```Horst```)
- absolute values (```centroid```)
- squared values (```factorial```)

Only, the horst scheme penalizes structural negative correlation. The factorial scheme discriminates more strongly the blocks than the centroid one.

### 2.5. Pre & post-analysis functionalities
Before running the RGCCA analysis, the optimal tau (for non-sparse models) and the sparsity coefficient (for sparse ones) can be evaluated by comparing the performance of RGCCA with random models: (i) with random permutation of the lines within each of the blocks to break the link between them; (ii) if supervised, with random sampling of the individuals used for modeling and validation of the prediction on the remaining ones (cross-validation). Several combinations of parameters (tau or sparsity coefficient) are evaluated. The parameters of the best model can be manually entered into the model.

After running the RGCCA analysis, a ```Run bootstrap``` button appears at the bottom of the analysis panel. If the analysis is supervised, the user could ```Evaluate the model``` (with 5 k-folds; see after).

These analysis are launched after clicking on their button. New graphical tabs will appear for each of these analysis. To directly visualize the outputs, see the [last section](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#4-visualize-the-plot-tabs).

#### 2.5.1. Number of bootstraps, permutations, folds 
By default, the number of bootstraps and permutations is set to 10 runs and could be increased up to 1000 runs.
In cross-validation, the evaluation of the performance of the model is based on two phases: the first consists of training the model with a first set of data, the second of evaluating it on a second set of data independent of the first. There are several ways of dividing its original data set into a training set and an evaluation set. One of these is the ```k-fold```. This consists in separating its data set into k equal parts. Only one part will be used for the evaluation phase and the other k-1 parts for the training phase. In the same way, each of the parts will be used for a new evaluation and the others for a training of the model. By default, a 5-fold is performed.


## 3. Graphical parameters ('Graphic' parameter tab)

This parameter tab is observed only with the ```samples```, ```corcircle``` and ```fingerprint``` plot tabs (**Fig. 3**).

![graphic](../../img/shiny/graphicPar.png)

*Fig. 3 : When the samples tab is selected a graphical option panel appears, that includes: (i) the possibility to hide/print the names of the variables, (ii) the selection of the block to visualize, (iii) the components used in the plots, (iv) the loading of groups of response to color the samples and (v) a button to save all the plot in the folder of the Shiny application. In this example, the agriculture block will be selected as the block for the Y-axis.*

### 3.1. Display names
If activated (by default), the ```display names``` option shows the name of the points in the ```samples``` and ```corcircle``` plot tabs. If disabled, shapes are shown instead of text: one per group of modality (see [section 3.4.](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#34-color-the-samples)).

### 3.2. Block (for the x/y-axis)
By default, if selected, the ```samples```, ```corcircle``` and ```fingerprint``` plot tabs are shown with the ```superblock``` (i.e., the concatenation of all blocs; see [section 2.3.2](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#232-superblocks)) to visualize all the blocs together. If this option is disabled, by default, the last blocks in the drop-down menu ```block``` is used. In this menu, choose another block (e.g., ```agriculture```) to update the plots with your selection (**Fig. 3**). For the ```samples``` tab only, a ```block for the x-axis``` and a ```block for the y-axis``` could be selected.

### 3.3. Components (for the x/y-axis)
The ```component``` of the analysis allows to choose the space where the points are visualised. For ```samples``` and ```corcircle``` biplots tabs, either ```component for the x-axis``` or ```component for the y-axis``` could be set. By default, they are respectively set to the first and the second components. Their choices are limited by the number of components selected in the analysis (defined in the 2.2. section). If the number of components in RGCCA were higher than two (not allowed in the Russet example, because of the industry block), the ```component for the x-axis```, for example, could be set to the third one.

### 3.4. Color the samples
On the ```samples``` plot tab only, one could select a variable to color the points according to a response. For this, load the ``` political_system.tsv``` file in the corresponding ```groups of modalities``` box to update the plot. he expected format is a CSV file tabulation-separated with: 
- qualitative or quantitative values (decimals separated by '.') with missing values labelled as "NA"
- samples in lines, labelled in the 1rst column with the same sample names as the blocks (some samples could be missing)
- a header containing the names of the columns

### 3.5. Number of top variables
Available only for ```fingerprint``` plot tabs, the ```number of top variables``` is automatically set to the number of variables in the selected blocks until a maximum value of 100. For example, with Russet data, eleven "top" variables could be visualised by default on the superblock.

### 3.6. Save the graphics
All graphics could be saved in the tool folder by running the ```Save all``` button (see the [description of the outputs](https://github.com/rgcca-factory/RGCCA/tree/master#output-files)).

### 3.7. Dynamic graph parameters header
In ```samples```, ```corcircle``` and ```fingerprint``` tabs, a header will appear allowing to dynamically explore the graph (see [**Fig. 6**](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#43-samples)). From left to right: 
- camera icon saves the plot
- magnifying icon zooms in the plot
- crossed arrows moves the plot
- home icons resets all the changes

### 3.8. Bootstrap indexes for x- and y- axis
These options allow to customize the choice of parameters to be displayed in the X or Y axis.
- ```bootstrap ratio```: average of the bootstrap weights divided by their standard error.
- ```significant 95% interval``` (disabled for SGCCA): weights where 95% of the values do not fluctuate around zero (no sign change).
- ```occurrences``` (available for SGCCA only): non-zero weight proportion
- ```mean boostrap weights```
- ```RGCCA weights```

### 3.9 Display cross-validation
If cross-validation is launched, this option shows or hides its predictions on the ```samples``` plot tab (see [**Fig. 6**](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#43-samples)).

## 4. Visualize the plot tabs

Please, make sure that these options are set to visualise the same plots than those in the next examples:
- RGCCA method
- two components
- block scaled
- tau = 1
- a superblock
- factorial scheme

After clicking on the ```run analysis``` button, a set of graphical tabs will appear on the right. By navigating between them, the user could visualize the output of the analysis in various format. On this plots, for the ```samples```, ```corcircle``` and ```fingerprint``` tabs, for each axis of the selected block, the corresponding percent of average explained variance is indicated. Two additional plot tabs could appear if additional analyses are carried out.

Every modified data or analysis parameter resets the plot tabs and the analysis needs to be re-run to apply the changes.

### 4.1. Connection between blocks
The first tab summarizes the connection between each block: a link corresponds to a "1" value, in the matrix connection file (**Fig. 4**; see [section 2.3.1.](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#231-loading-a-connection-file)). For each block:
- "P" is the number of variables
- "N" is the number of lines (here, each block has the same number of line)
- "tau" is the shrinkage parameter and "sparsity" is the sparsity coefficient (see the [2.4.1 & 2.4.2 sections](https://github.com/rgcca-factory/RGCCA/blob/master/inst/shiny/tutorialShiny.md#24-other-rsgcca-parameters)). The tau parameter could be shown for each component if the optimal option is selected

![connection](../../img/shiny/connection.png)

*Fig. 4 : Connection between each block of the RGCCA and the superblock with 47 common rows between blocks*

### 4.2. Average variance explained (AVE)
In the second tab on the right panel, the average variance explained (AVE; in X-axis) is represented in percent for each block (in Y-axis) and each component (one color per component) (**Fig. 5**). The subtitle informs about the AVE for the two first of the outer model (weighted average of the AVE of each block).

![ave](../../img/shiny/ave.png)

*Fig. 5 : Average variance variance explained (in %) for each block and for the two first components of the RGCCA*

### 4.3. Samples
The first tab is the projection of the sample coordinates in the selected component of the analysis and, by default, on the
superblock (a concatenation of all the blocks) (**Fig. 6**). If a ```response``` file is loaded, each sample is colored according to this variable. In the Russet example, the X-axis could discriminate a dictatorship (with upper values on this axis than the two other political systems), whereas the Y axis discriminates an unstable democracy (with upper values than the others).

![individuals](../../img/shiny/samples.png)

*Fig. 6 : Samples coordinates on the two first components for the superblock of the RGCCA after loading the " political_system.tsv" file. By selecting the option in the header of the plot, a zoom could be done on a selected part of the graph (e.g., the upper right part). Then, all the modifications could be reset with the "house" icon.*

When a cross-validation has been carried out and the graph option has been activated, the graph of individuals allows the actual values of the GCCA (in red) to be contrasted with those predicted by the validation (in green). 
Here, by selecting a supervised mode with industry as the default response block, we can see that, except for Iraq, the prediction of the leave-one-out is close to reality. The prediction score is the average of the absolute values of the differences between the actual and predicted data centered and reduced. At the start of the cross-validation in the analysis tab, a pop-up indicates a value of 0.0919. This figure can be seen as a percentage difference between predicted and actual data.

![cross-validation](../../img/shiny/cross-validation.png)

*Fig. 7 : With a cross-validation, it is possible to compare his predictions (in green) with the real coordinates of the individuals (in red).*

### 4.4. Corcircle
The second one corresponds to the Pearson correlation between the variables of the block and the selected components in the analysis (by default, on the two first components) (**Fig. 7**). The circle is a 1 correlation and the dotted one is a 0.5 correlation. If the superblock is selected, colors correspond to the belonging of each variable to each block. Only the 100th variables the most correlated to each axis are printed.

![corcircle](../../img/shiny/corcircle.png)

*Fig. 8 : Correlation between each variable of each block (by using the superblock) and the two first components of the RGCCA. The "print names" parameter in the graphical parameter panel is disabled. By flying over a point, an informational label appears. This label gives the name and the X and Y coordinates of the point.*

### 4.5. Top variables
The next tab also represents the same correlation of the variable with the selected component (on the X-axis; 1 by default). The top variables are sorted decreasingly (on the Y-axis) in a histogram among the selected block (superblock, by default) (**Fig. 8**). 

![top_variables](../../img/shiny/fingerprint.png)

*Fig. 9 : Top 11 variables among all the blocks (by using the superblock) with higher correlation with the first component of the RGCCA. "Gnpr" (belonging to the industry block) shows a correlation of 0.859 with this component.*

Here, "labo" from the industry block (the amount of labor force in agriculture) is the variable the most positively correlated to the X-axis. On the opposite, "gnpr" (gross national product) from the industry block and "demostab" (stable democracy) from the political block are the most negatively correlated variables. In other terms, these variables are the most importants on the first component of the RGCCA. Countries with an unstable democracy are more associated with a lower "rent" (percent of farmers that rent their land). Otherwise, those with a dictatorship system are more associated with a higher labor force in agriculture values and a less gross national product (and inversely for the stable democracy case).

### 4.5. Bootstrap

The graph represents the bootstrap ratio (average of the bootstrap weights divided by their standard error) on the X-axis according to their significance (*) or not (ns). 
On this graph, the gross national product is the variable with the highest average weight among the bootstraps. Having a negative weight in the original RGCCA, it is coloured red.
The unstable democracy and the land rent fluctuate around 0 and are therefore not considered positive. They have on average a very low and positive weight in the original RGCCA (in green).

For this type of analysis and the next (permutation), due to a random process, different results can be obtained between two runs. It is advisable to progressively increase the number of runs until the results converge. Here, 100 bootstraps could be enough (~20 seconds to run 1000 bootstrap).

![bootstrap](../../img/shiny/bootstrap.png)

*Fig. 10 : The highest average weight, crossed with the significance, allows the selection of the gross national product as a variable of interest.*

### 4.6. Permutation
Here, penalties are selected as the choice of parameters to permute and a number of permutations of 1000 (could be more than 5 minutes to run). Be careful with large data sets: test the time it takes 10 permutations before increasing the number of permutations. On this graphics, the abscissa shows the combinations tested during the permutations. On the ordinates are the Z-score or significance score. The three horizontal bars in dashed lines correspond to the three significance thresholds from the lowest to the highest bar: 95%, 99% and 99.9%.

In our example, with a risk of 99.9%, the best choice of combinations is the one with the lowest sparsity coefficients. This choice can guide our decision in setting the parameters of our model.

![permutation](../../img/shiny/permutation.png)

*Fig. 11 : A permutation allows to select a set of significant parameters, here the combination of the lowest sparsity coefficients.*
