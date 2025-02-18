# BackupGene

## Overview

Bioinformatic pipeline designed for gene compensatory (or backup genes) prioritization, supporting various execution modes for data preparation, ranking, integration, and report generation. The mode of execution is controlled by the `exec_mode` parameter, determining the stage of the process.

## Execution Modes

The workflow is controlled from `daemon.sh` and provides several execution modes (`exec_mode`) which determine the actions performed. Below is a description of each mode:

### `download_translators`
- Downloads translator tables for gene annotations.
- Creates a `./translators` directory and fetches files like `HGNC_symbol` from external resources.

### `control_preparation`
- Prepares control gene data for further analysis.
- Handles two types of control preparation based on the second argument (`bk` or `targets`):
  - `bk`: Prepares control genes using backup data.
  - `targets`: Aggregates control genes using Open Targets data.

### `control_type`
- Defines control gene types, such as robust and non-robust genes.

### `ranking`
- Computes rankings based on non-integrated similarity kernels.
- Executes ranking for multiple kernels and annotations.
- Saves results in the `rankings` directory.

### `integrated_ranking`
- Similar to `ranking`, but involves integrated kernels.
- Uses different integration types (mean, median, max) and stores results in the `integrated_rankings` folder.

### `report`
- Generates a report based on computed results.
- Creates an HTML report using ranking data, CDFs, and additional control gene information.
- Requires Python environments and templates for report creation.

### `check`
- Validates the execution flow in a specific folder.
- Logs the status of all executions in that folder.

### `recover`
- Recovers data from a prior execution.
- Reprocesses data in case of failed executions.

## Input Variables

- `exec_mode`: Determines the mode of execution.
- `add_opt`: Additional options specific to the execution mode.
- Various paths are configured, including input directory, output folders, and scripts directories. It is worh to mention that, as this reposiotory depends from results extracted in [GraphPrioritizer](https://github.com/federedef/GraphPrioritizer.git), in `output_folder_GraphPrioritizer` variable, you can define the path to these results. 

## Path Setup

- `input_path`: Current working directory.
- `PATH`: Includes directories for necessary scripts and programs.
- `autoflow_scripts`: Path for AutoFlow-related scripts.
- `daemon_scripts`: Path for daemon scripts.
- `control_genes_folder`: Folder containing control gene data.
- `output_folder`: Main directory for storing results.

## Custom Variables

- `annotations`: List of annotations, such as disease, phenotype, and biological processes.
- `kernels`: List of kernels for similarity calculations.
- `integration_types`: Types of integration strategies (mean, median, max).

## Functions and Workflow Stages

The script utilizes `AutoFlow` to execute certain stages of the workflow, using specific arguments for each task.

## Output

The output is organized in multiple directories such as `rankings`, `integrated_rankings`, and `report`. The generated reports are stored in HTML format and accompanied by images and data files.

## Dependencies

- `AutoFlow`: Used for executing the ranking and integration stages.
- `wget`, `awk`, `grep`: Command-line utilities for data processing.
- Python: Required for report generation and other processing steps.

## Notes

- As mentioned before, this repository is dependent from results obtained with [GraphPrioritizer](https://github.com/federedef/GraphPrioritizer.git).
- The script supports various stages that can be customized based on the dataset and the desired analysis.

## Example Usage

CodeIn
./script.sh ranking
./script.sh control_preparation targets
./script.sh report final_report_name
CodeOut

This will execute the respective stages based on the provided `exec_mode`.
