# QMSFC

Data and code for the `Multiomics surface receptor profiling of the NCI-60 tumor cell panel uncovers novel theranostics for cancer immunotherapy` paper.

## Clone this repository

We can clone via:

```bash
git clone git@github.com:qbicsoftware/QMSFC.git
```

Then change your directory:

```bash
cd QMSFC
```

## Install basic required software

Some basic software is required in order to run the R analysis scripts. Please install `Rstudio` as described in https://www.rstudio.com/products/rstudio/download/#download. To get a `miniconda` installation please follow https://conda.io/projects/conda/en/latest/user-guide/install/index.html.

## Prepare conda environments QMSFC and QMSFC_2

Once you have installed `miniconda` you can add the QMSFC environment as follows:

```bash
conda env create -f QMSFC_environment_21042022.yml
```

For the QMSFC_2 environment, you have to execute:

```bash
conda env create -f QMSFC_environment_22042022.yml
```

## Run quality control analysis

Ensure that you are in the QMSFC environment:

```bash
conda activate QMSFC
```

For `Rstudio` to work properly, we tell it the path to th R version used in the conda environment:

```bash
export RSTUDIO_WHICH_R=~/miniconda3/envs/QMSFC/bin/R
```

If you installed `miniconda` to a different location than your home folder, please adjust the path above.

Now open `Rstudio` from the command line:

```bash
rstudio
```

Depending on your OS, it may also be `rstudio-bin`.

Navigate to `File` and click on `Open File...`. Navigate to the `IsotypeControl_QC` folder and open `IsotypeControl_QC.R`. Select all lines and press `SHIFT+ENTER` in order to re-run the analysis.

## Run pairwise correlation analysis for cell lines measured in 2 weeks

Ensure that you have the QMSFC environment activated and `Rstudio` is open as in [Run quality control analysis](#run-quality-control-analysis).

Navigate to `File` and click on `Open File...`. Navigate to the `week1_week2` folder and open `week1_week2.R`. Select all lines and press `SHIFT+ENTER` in order to re-run the analysis.

## Run homogenize identifiers script

Ensure that you have the QMSFC environment activated and `Rstudio` is open as in [Run quality control analysis](#run-quality-control-analysis).

Navigate to `File` and click on `Open File...`. Navigate to the `Homogenize_Identifiers` folder and open `harmonize_IDs.R`. Select all lines and press `SHIFT+ENTER` in order to re-run the analysis.

## Run Tx curation

Ensure that you are in the QMSFC_2 environment:

```bash
conda activate QMSFC_2
```

For `Rstudio` to work properly, we tell it the path to th R version used in the conda environment:

```bash
export RSTUDIO_WHICH_R=~/miniconda3/envs/QMSFC_2/bin/R
```

If you installed `miniconda` to a different location than your home folder, please adjust the path above.

Now open `Rstudio` from the command line:

```bash
rstudio
```

Depending on your OS, it may also be `rstudio-bin`.

Navigate to `File` and click on `Open File...`. Navigate to the `Curate_ArrayData_April2022` folder and open `curate_array_April2022.R`. Select all lines and press `SHIFT+ENTER` in order to re-run the analysis.

## Run Px curation

Ensure that you have the QMSFC_2 environment activated and `Rstudio` is open as in [Run Tx curation](#run-tx-curation).

Navigate to `File` and click on `Open File...`. Navigate to the `Curate_Px_April2022` folder and open `curate_px_gholami_April2022.R`. Select all lines and press `SHIFT+ENTER` in order to re-run the analysis.

## Run MCIA analysis

Ensure that you have the QMSFC environment activated and `Rstudio` is open as in [Run quality control analysis](#run-quality-control-analysis).

Navigate to `File` and click on `Open File...`. Navigate to the `multiomics_March2022_NO_LOX_NO_ISOTYPE_CONTORLS` folder and open `omicade4_March2022_NO_LOX_NO_ISOTYPE_CONTROLS.R`. Select all lines and press `SHIFT+ENTER` in order to re-run the analysis.

## Run Recount2 analysis

For each of the folders starting with `recount2` you can run the same analysis steps. In the following an example for COLON.

Ensure that you have the QMSFC_2 environment activated and `Rstudio` is open as in [Run Tx curation](#run-tx-curation).

Navigate to `File` and click on `Open File...`. Navigate to the `recount2_co_April2022` folder and open `curate_px_gholami_April2022.R`. Select all lines and press `SHIFT+ENTER` in order to re-run the analysis.
