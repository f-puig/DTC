This Markdown document presents a pipeline for the investigation of 1H NMR spectra with the *Decision Tree of Correlation - Multivariate Curve Resolution - Alternating Least Squares* (or **DTC-MCR-ALS**) method.

The **DTC-MCR-ALS** method is a powerful tool for deciphering the composition of metabolomics samples from 1H NMR spectra, since it can resolve the pure NMR spectra and the concentrations of the metabolites found in these metabolomic samples.

The **DTC-MCR-ALS** method results very useful when the metabolic composition of the samples to be analyzed is not known in advance. Since resonances are investigated without having been assigned first, this method can be considered as an **untargeted** method for the analysis of 1H NMR metabolomics samples.

With this method, the overlapped resonances are separated using the power of [MCR-ALS](https://mcrals.wordpress.com/). **MCR-ALS** is a **curve resolution** method, and it does not rely on [deconvolution](https://en.wikipedia.org/wiki/Deconvolution) nor [curve fitting](https://en.wikipedia.org/wiki/Curve_fitting). For this reason, the best results are obtained for cases with **high metabolic variance** among the studied samples, even if resonances present  **high degree of overlap**. More information regarding the **DTC-MCR-ALS** method can be found in this **[article](https://doi.org/10.1016/j.aca.2017.02.010)**.

The **DTC-MCR-ALS** method has been implemented under **MATLAB** environment as a set of functions, which can be downloaded from this respository and from the [MCR-ALS](https://mcrals.wordpress.com/) webpage.

A **guide** for the execution of the DTC-MCR-ALS method is given below. 



The steps in the execution of the DTC-MCR-ALS method can be summarized in the following four:

**1)** Windowing of the 1H NMR dataset.

**2)** MCR-ALS resolution of the NMR spectral windows.

**3)** Grouping the resolved profiles in **step 2** using the DTC method.

**4)** MCR-ALS resolution of the whole 1H NMR dataset.



## 0. 1H NMR DATASET USED AS EXAMPLE ##

In this **guide**, a 1H NMR metabolomics dataset is used as an example. This dataset can be found in the **NMRdataset.mat** file. 

The given dataset consists of 60 simulated 1H NMR spectra of mixtures composed by 10 metabolites. In these simulated spectra, random noise was added to simulate the background signal observed tipically in the real 1H NMR spectra. The used concentration values were obtained from a previous metabolomic experiment using real NMR samples of yeast extracts (see **[article](https://doi.org/10.1016/j.aca.2017.02.010)**).

In **Figure 1** below, selected regions from these samples are shown.

![Fig2](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig2.png)

***Figure 1***. *Highlighted regions from the simulated 1H NMR metabolomics dataset.*



In the following plot in **Figure 2**, relative concentrations used for the set of 60 samples.

![Fig1](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig1.png)

***Figure 2***. *Relative concentration values for the 10 metabolites used to simulate the 1H NMR metabolomics dataset.*



More information about this dataset can be found in this **[article](https://doi.org/10.1016/j.aca.2017.02.010)**. In the article, this dataset is referenced as dataset **X2**.

 

## 1. Windowing of the 1H NMR dataset ##

The first step of the DTC strategy consist in the reduction of the complexity of the dataset. This simplification is achieved by analyzing the dataset in small pieces instead of as a whole.

Thus, the dataset is split in smaller datasets, each one covering a small spectral window. The number of spilts must be the highest possible, but avoiding to split a resonance between windows. For this reason, the limits of these windows must be located in chemical shift positions that are only descriptive of background noise.

If two overlapped regions are slightly collided (see **Figure 3**), then it is preferable to set the limit for these windows in the chemical shift position associated with the minimal intensity (blue arrow in **Figure 3**).

![Fig3](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig3.png)

***Figure 3.*** *Defining spectral windows limits in overlapped regions.* Two distinct overlapped regions (**1** and **2**) can be observed in the figure. The blue arrow points to the chemical shift position that will be used to separate the two overlapped regions.

Therefore, depending on the degree of overlapping, the spectral data contained in each window will be:

- An isolated resonance.
- A set of overlapped resonances.

For this small example dataset, only 26 windows were needed, but for datasets that contain more resonances, more windows will be needed. For instance, for a 1H NMR dataset of orine samples, which can contain more than 500 resonances from hundreds of compounds, more than 150 windows may be needed.

In **Figure 4** below, the suggested limits for the analysis of the example dataset are shown.

![Fig4](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig4.png)

***Figure 4.*** *Limits used for the analysis of the example dataset.*

The list of limits suggested for this example are provided in the ``limits`` matrix in the file **limits.mat**. In this matrix, for each window, the starting and ending positions are given. For example, window 1 (**w_1**) comprises the set of intensities between the **1st** position (ẟH= 0.8 ppm) and the **547th** (ẟH= 1.2 ppm). 

After having stablished these limits, the dataset must be split into spectral windows. This can be simply achieved with the function **splitNMR.m**.

```matlab
windows = splitNMR(NMR,limits);
```

This function creates a **structure** named ``windows`` that contains the set of spectral windows generated from the original dataset.



## 2. MCR-ALS ##
#### 2.1. MCR-ALS resolution of the NMR spectral windows

Each one of the generated spectral windows must be now resolved with MCR-ALS.  The MCR-ALS chemometric method is able to separate in different components those resonances that show different variance among samples.

From every MCR-ALS analysis, a matrix with the resolved concentration profiles, ```copt```, and a matrix with the resolved spectral profiles, ```sopt```, are generated. 

MCR-ALS analysis of every window can be performed with the following command:

```matlab
[copt_w_i,sopt_w_i] = als_windows(windows_i,num_comp);
```

In this MCR-ALS, the number of components must be set by the user. Thus, it is convenient to repeat the MCR-ALS analysis with a different number of components in order to determine the optimal number.

For instance, given the **w_10** spectral window in **Figure 5**, the MCR-ALS resolution using 1 and 2 components is shown in **Figure 6** and **Figure 7**, respectively. We can appreciate that the optimal number of components is 2 because this spectral window is composed of 2 triplets.

![Fig5](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig5.png)

***Figure 5.*** *Windows 10.*

![Fig6](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig6.png)

***Figure 6.*** *MCR-ALS resolution of windows 10 using 1 component. In the left, the resolved C profile. In the right, the resolved ST profile.*

![Fig7](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig7.png)

***Figure 7.*** *MCR-ALS resolution of windows 10 using 2 componenta. In the left, the resolved C profiles. In the right, the resolved ST profiles.*

For regions where peaks presents a low signal-to-noise ratio, it is advisable to add additional components to explain noise in order to resolve better defined resonances.

For the given example dataset, **33 components** were needed to resolve all the different spectral windows that comprise the whole NMR dataset.



#### 2.2. Selection of the rellevant components

Components that explain noise must be removed. To determine which components should be excluded, the **evaluate_features.m** Matlab function can be employed:

```matlab
[ignored] = evaluate_features('file.mat');
```

This function reads a ``.mat`` file (here, named ``file.mat``) containing the set of ```copt``` and ```sopt``` matrices obtained in the previous step and plots them (see **Figure 8**). In these plots, components are labelled according to the total number of components resolved in all the windows. For instance, for **windows 5**, 2 components were used. These two components are, according to this labelling, component number 7 and component number 8, respectively (see **Figure 8**).

![evaluate_feature](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\evaluate_feature.png)

***Figure 8.*** *Resolved components for windows 5. In the left, the resolved C profiles. In the right, the resolved ST profiles.*

To create the ```.mat``` file that will be read by the **evaluate_features.m** function, variables (```copt```and ```sopt```) must be selected from the workspace, and saved together after pressing the right button with mouse and click '*save as*'. To distinguish the different ```copt```and ```sopt``` variables obtained from the different windows, the suffix 'w#' (# is the number of the window) will be used. For instance, the resolved spectra profiles for **w15** are saved in the variable 'sopt_**w15**'.

The **evaluate_features.m** Matlab function also checks that all variables used in the DTC (see **step 3**) have the expected size. For this reason, the  ``.mat`` must also contain the vector with the ppm values (which must be named ```ppm```).

In the given example, after application of the **evaluate_features.m** function, 2 components were considered to be not descriptive of resonances. These components were component 7 and component 21. The output of this function contains the indeces of the components that will not be used in the next steps of the analysis.

```matlab
[ignored] = [7 21];
```

Thus, after ignoring these 2 components, only **31 components** (or features) were used in the DTC analysis.



#### 2.3. Preparing the data for the DTC

Now that the set of rellevant components are known, the next step is to unify all the matrices descriptive of **C** and **ST** profiles into just two variables (one for the set of **C**, and one for the set of **ST**) in the adequate format for the function responsible of the DTC (**selNMRfeat.m**). This operation is performed using the ***supersopt_supercopt_ppm.m*** Matlab function.

```matlab
[supersopt,supercopt,names_sopts] = supersopt_supercopt_ppm('file.mat', ignored);
```

The inputs of this function are the ``'file.mat'``  used in the **evaluate_features.m** function (see **step 2.2**), and the ``ignored`` vector obtained as output of the same function. In the case that all the components were rellevants (0 components describing noise), then only the ``'file.mat'`` filename must be given.

```matlab
[supersopt,supercopt,names_sopts] = supersopt_supercopt_ppm('file.mat');
```

The outputs of this function are:

- ``supersopt``: a cell array containing all the rellevant **ST** profiles. *supersopt* has as many cells as number of rellevant components. 
- ``supercopt``: an augmented matrix containing all rellevant **C** profiles. The dimensions of *supercopt* are *m* by *n*, being *m* the number of samples, and *n* the total number of rellevant components. 
- ``names_sopt``: the names of the variables used to build up the *supersopt* cell array.



##  3. Grouping the resolved profiles using the DTC method ##

The different features corresponding to the same metabolite are grouped using the DTC method. With this method, these features are grouped one by one, from the most correlated ones to the least correlated ones. This grouping strategy is human-driven, and it can be easily performed using the **selNMRfeat.m** Matlab function:

```matlab
[group_indeces]=selNMRfeat(NMR, ref, ppm, supercopt, supersopt, grouped_indeces, new);
```

The inputs for this funtion are the following:

- ``NMR``: the original 1H NMR dataset. 
- ``ref``: a reference 1H NMR spectrum.
- ``ppm``: the vector of chemical shift positions.
- ``supercopt``: the *supercopt* matrix generated in **step 2.3**.
- ``supersopt``: the *supersopt* cell array generated in **step 2.3**.
- ``group_indeces``: the cell array containing the indeces of the grouped features.
- ``new``: reflects whether the next feature to be picked must be grouped with an existent group (``new = 0``), or a new group of features is started (``new = 1``).



Before running the **selNMRfeat.m** function, the ``group_indeces`` cell array needs to be initiallized first:

```matlab
group_indeces = {};
```



And the DTC can be finally started:

```matlab
[group_indeces]=selNMRfeat(NMR, 1, ppm, supercopt, supersopt, grouped_indeces, 1);
```

After introducing this command, a 4-plot window will be opened (**Figure 9**), and the following message will be displayed:

> Higher correlation is found between 6 and 13 
> Correlation value is 0.999963 
> Group NMR feature: 1-YES/0-NO

![DTC1](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\DTC1.png)

**Figure 9**.  DTC output.

The 4 plots generated in this step contain different data that can guide the user in the decision of whether the two features belong to the same compound or not:

- The first plot (upper): a plot with all 1H NMR spectra overlaid.
- The second plot: a plot showing the selected features for the group being built. The resonance in red is the one being considered to be grouped or not with the previously grouped features (in blue).
- The third plot: a plot with the spectral features shown in the second plot, scaled by the concentrantrations estimated in all the samples. Every different color represents one sample.
- The fourth plot (bottom): for every sample, the ratio between the concentration associated to the feature to be selected and one of the concentrations of the already grouped features. If this value is steady for all samples, then it is implied that the two resonances are relative to the same compound



Since the two selected features present a high degree of correlation (0.999963), their variance is similar among the analyzed samples (third plot), and the concentration ratio between these two features is mainly maintained for all the samples, we can conclude that **they are from the same metabolite** and they must be grouped toguether. In order to group them, we type **1** on the command line:

> Group NMR feature: 1-YES/0-NO **1**



After this, the **selNMRfeat.m** function will end, and the ``group_indeces`` cell array will contain in the first position of the cell the indeces of the two grouped features, 6 and 13.

Now, the DTC process must be continued (until finishing the grouping of all features). To continue, we must write the same line as before, but using ``new=0``.

```matlab
[group_indeces]=selNMRfeat(NMR, 1, ppm, supercopt, supersopt, grouped_indeces, 0);
```

With this command, the feature to be picked will be the one that shows the highest correlation to either feature 6 or to feature 13.

These are the outputs:

> Higher correlation is found between 6 and 2 
> Correlation value is 0.999957 
> Group NMR feature: 1-YES/0-NO

![DTC2](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\DTC2.png)

**Figure 10**.  DTC output.

From this second representation, we can guess that feature 2 is also part from the same group of features as the featurs 6 and 13, and so we can add it with the other two.

> Group NMR feature: 1-YES/0-NO **1**

The process is continued, using the very exact command line:

```matlab
[group_indeces]=selNMRfeat(NMR, 1, ppm, supercopt, supersopt, grouped_indeces, 0);
```

> Higher correlation is found between 6 and 4
> Correlation value is 0.999431
> Group NMR feature: 1-YES/0-NO

![DTC3](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\DTC3.png)

**Figure 11**.  DTC output.

From the representation in **Figure 11**, we can conclude that the feature 4 also is associated to the group of features being built. 

> Group NMR feature: 1-YES/0-NO **1**

The process is continued, using the very exact command line:

```matlab
[group_indeces]=selNMRfeat(NMR, 1, ppm, supercopt, supersopt, grouped_indeces, 0);
```

> Higher correlation is found between 2 and 18 
> Correlation value is 0.9888451
> Group NMR feature: 1-YES/0-NO

![DTC4](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\DTC4.png)

**Figure 12**.  DTC output.

From the representation in **Figure 12**, we can conclude that the feature 18 also is associated to the group of features being built. 

> Group NMR feature: 1-YES/0-NO **1**



If this process is continued, the new suggested pair is:

> Higher correlation is found between 2 and 18
> Correlation value is 0.946835
> Group NMR feature: 1-YES/0-NO

![DTC6](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\DTC6.png)

**Figure 13**.  DTC output.

By looking carefully to the inputs (**Figure 13**), we can appreciate that the spectral profile of the picked resonance (multiplet centered at 2.1 ppm) do not follow the same pattern for all the samples (third plot in **Figure 13**). This can be detected because the sample with the highest feature is green, whereas this is not the case for the other three features shown in the same plot. In fact, for threse three other features (multiplet at 1.5, 1.7, and 1.9 ppm), samples with resonances in green are much smaller than the others.



Therefore, since this feature should not be grouped with the group being built, we must not add anymore feature to the group being build up:

> Group NMR feature: 1-YES/0-NO **0**



In this point, ```group_indeces``` is a cell object with 2 entries:

- ```group_indeces{1,1}= [ 2 4 6 13 18 ];```
- ```group_indeces{1,2}= [ 7 ];```

However, since feature ``7`` was not grouped in the first group, ``group_indeces{1,2}`` must be deleteted in order to be considered to be part of another group. This can be directly performed by deleting the second column of the cell array:

![delete_columns](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\delete_columns.png)

**Figure 14**.  Deleting columns.



After having removed ``group_indeces{1,2}``, we can now proceed to continue with the DTC with the remaining features. Since a new group of features need to be started, we need to use ```new=1```:

```matlab
[group_indeces]=selNMRfeat(NMR, 1, ppm, supercopt, supersopt, grouped_indeces, 1);
```

> Higher correlation is found between 27 and 29 
> Correlation value is 0.998956
> Group NMR feature: 1-YES/0-NO

As it can be seen, the chosen features show a much higher correlation (**0.998956**) than the last features picked in the former group of features.![DTC5](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\DTC5.png)

**Figure 15**.  DTC output.

And this process is being performed until all the features have been grouped. With the given example dataset, features should be grouped in the following groups (metabolite name given inside the parenthesis):

- Group 1 (L-lysine): 2, 4, 6, 13, 18.
- Group 2 (L-histidine): 14, 15, 16, 21, 27, 29.
- Group 3 (L-glutamic acid): 7, 8.
- Group 4 (L-ornithine): 3, 5, 12.
- Group 5 (AMP): 23, 24, 25, 30, 31.
- Group 6 (L-asparagine): 10, 11, 22.
- Group 7 (L-tyrosine): 20, 26, 28.
- Group 8 (Citric acid): 9.
- Group 9 (Glycine): 17.
- Group 10 (L-leucine): 1.
- Group 11 (mixture): 19.

From the DTC strategy, 11 groups of resonances were obtained, but only 10 of them can be associated to real metabolites, while the other group contains resonances that could not be properly grouped because they were poorly resolved in the MCR-ALS due to amiguity and rank defiency problems. In other words, resonances contained in the analyzed windows had similar variance among the studied samples, affecting negatively in the resolution. Nevertheless, it should be noted that all the metabolites could be (at least partially) identified by using any of the other 10 group of metabolites.

The next step is to perform a MCR-ALS analysis on the set of selected resonances, which is covered in the following section below.



## 4. MCR-ALS resolution of the whole 1H NMR dataset

This last step is performed in order to obtain MCR-ALS components containing the spectral information of the different groups unified. More over, the unification process of the spectral features will produce, as a counterpart, that each component will be associated to only one concentration value per sample (instead of presenting one concentration value per sample and feature).

To carry out the MCR-ALS analysis, we need first to generate some inputs. These inputs are the matrices with the **initial estimates** and the one with the **spectral selectivity constraints**. These inputs are automatically generated after application of the **nmr_ssel_iniest.m** Matlab function.

```matlab
[ssel, ini_est]=nmr_ssel_iniest(supersopt, supercopt, grouped_copt_indeces, ppm);
```



Finally, the only last thing to do is to call the MCR-ALS method itself:

```matlab
[r2,copt,sopt]=als_command_ssel_dtc(d,ini_est,nit,tolsigma,ssel);
```

The inputs for this MCR-ALS method are:

- ``d``: the experimental data matrix.
- ``ini_est``: the matrix with the initial estimates.
- ``nit``: the number of iterations.
- ``tolsigma``: the convergence criterion.
- ``ssel``: the matrix with the spectral selectivity constraints.

For this case, we used 1000 as the number of iterations and 0.1 as the convergence criterion:

```matlab
[r2,copt,sopt]=als_command_ssel_dtc(NMR,ini_est,1000,0.1,ssel);
```



Up to a certain point, the iterative MCR-ALS method will reach the convergence state and it stops. With this, a matrix containing the 1H NMR spectra (``sopt``, covering the full spectral range) and another containing the relative concentrations for the corresponding metabolites (``copt``) will be obtained. 

![Fig16](C:\Users\putxv\Desktop\github\DTC-MCR-ALS\Fig16.png)

**Figure 16**.  MCR-ALS analysis.

## References: ##

- Puig-Castellví F, Tauler R, Alfonso I (2017). Untargeted assignment and automatic integration of 1H NMR metabolomic datasets using a multivariate curve resolution approach. Analytica Chimica Acta. https://doi.org/10.1016/j.aca.2017.02.010.

## Contact: ##

puig.francesc@gmail.com
