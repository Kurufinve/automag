# Magnetocrystalline Anisotropy Energy Calculation

<cite>
**Referenced Files in This Document**   
- [1_submit.py](file://4_mae/1_submit.py)
- [2_plot_results.py](file://4_mae/2_plot_results.py)
- [3_submit.py](file://4_mae/3_submit.py)
- [4_plot_results.py](file://4_mae/4_plot_results.py)
- [MAE.py](file://4_mae/MAE.py)
- [validate_consistency.py](file://4_mae/validate_consistency.py)
- [input_example_enhanced.py](file://4_mae/input_example_enhanced.py)
- [IMPLEMENTATION_SUMMARY.md](file://4_mae/IMPLEMENTATION_SUMMARY.md)
- [run_vasp.py](file://ase/run_vasp.py)
</cite>

## Table of Contents
1. [Introduction](#introduction)
2. [Two-Stage Workflow](#two-stage-workflow)
3. [Core Calculation Logic](#core-calculation-logic)
4. [Result Validation and Consistency](#result-validation-and-consistency)
5. [Input Configuration and Examples](#input-configuration-and-examples)
6. [Integration with FireWorks and VASP](#integration-with-fireworks-and-vasp)
7. [Troubleshooting Common Issues](#troubleshooting-common-issues)

## Introduction
The Magnetocrystalline Anisotropy Energy (MAE) calculation module in automag-1 is designed to determine the magnetocrystalline anisotropy of magnetic materials by comparing the energies of non-collinear density functional theory (DFT) calculations with different magnetization directions. This module plays a crucial role in understanding the magnetic properties of materials, particularly in identifying the easy and hard magnetization axes, which are essential for applications in magnetic storage and spintronics. The MAE is calculated as the energy difference between the most stable (lowest energy) and least stable (highest energy) magnetization directions, providing insights into the material's magnetic anisotropy. The module is structured to ensure robust and accurate results through a two-stage workflow, which includes initial convergence and validation followed by directional energy comparisons. This approach not only enhances the reliability of the calculations but also facilitates the identification of any inconsistencies or errors in the results.

**Section sources**
- [README.md](file://README.md#L1-L277)
- [MAE.py](file://4_mae/MAE.py#L1-L1086)

## Two-Stage Workflow
The MAE calculation module in automag-1 employs a two-stage workflow to ensure accurate and reliable results. The first stage involves initial convergence and validation, which is handled by the scripts `1_submit.py` and `2_plot_results.py`. These scripts are responsible for setting up and executing the initial non-collinear DFT calculations to determine the optimal k-point mesh and energy cutoff (ENCUT) values. The `1_submit.py` script submits the necessary VASP jobs to the remote database, while `2_plot_results.py` processes the output to generate plots of energy versus k-point mesh and ENCUT values, helping to identify the converged parameters. This stage is crucial for ensuring that subsequent calculations are performed with accurate and consistent parameters.

The second stage focuses on the directional energy comparisons, which are managed by the scripts `3_submit.py` and `4_plot_results.py`. In this stage, the module performs non-collinear DFT calculations for various magnetization directions, defined by the angles θ (theta) and φ (phi). The `3_submit.py` script submits the VASP jobs for these directional calculations, while `4_plot_results.py` analyzes the results to determine the magnetocrystalline anisotropy energy. This stage involves comparing the energies of different magnetization directions to identify the easy and hard axes of magnetization. The two-stage workflow ensures that the calculations are both efficient and accurate, providing a comprehensive understanding of the material's magnetic anisotropy.

**Section sources**
- [1_submit.py](file://4_mae/1_submit.py#L1-L216)
- [2_plot_results.py](file://4_mae/2_plot_results.py#L1-L500)
- [3_submit.py](file://4_mae/3_submit.py#L1-L241)
- [4_plot_results.py](file://4_mae/4_plot_results.py#L1-L263)

## Core Calculation Logic
The core calculation logic of the MAE module is implemented in the `MAE.py` script, which orchestrates the non-collinear DFT calculations and energy comparisons. The script begins by defining the grid points for the angles θ and φ, which are used to specify the magnetization directions. The number of grid points for each angle is configurable, allowing for a balance between computational efficiency and accuracy. The script then sets up the VASP parameters, including the energy cutoff (ENCUT), k-point mesh, and other relevant settings, which are critical for the convergence of the calculations.

The `MAE.py` script uses the `pymatgen` library to handle the crystal structure and magnetic moments, ensuring that the calculations are performed on the correct atomic configuration. The magnetic moments are initialized based on the specified values for each magnetic atom, and the script ensures that the non-collinear nature of the calculations is properly accounted for. The script then submits the VASP jobs for each magnetization direction, using the `sbatch` command to manage the job submission on a high-performance computing cluster. After the calculations are completed, the script reads the output files to extract the total energies and magnetic moments, which are used to compute the MAE.

The energy comparisons are performed by fitting the energy data to a functional form that describes the anisotropy energy as a function of the magnetization direction. The script uses the `scipy.optimize.curve_fit` function to fit the data to a model that includes terms for the first and second-order anisotropy constants (K1 and K2). The fitted parameters are then used to determine the easy and hard axes of magnetization, providing a quantitative measure of the material's magnetic anisotropy. The `MAE.py` script also generates plots of the energy versus the magnetization direction, which are useful for visualizing the anisotropy and identifying any inconsistencies in the results.

**Section sources**
- [MAE.py](file://4_mae/MAE.py#L1-L1086)
- [2_plot_results.py](file://4_mae/2_plot_results.py#L1-L500)
- [4_plot_results.py](file://4_mae/4_plot_results.py#L1-L263)

## Result Validation and Consistency
To ensure the reliability and consistency of the MAE calculations, the `validate_consistency.py` script is used to compare the results of the original `MAE.py` script with those of the ported modular implementation. This script helps identify and debug any discrepancies between the two implementations, ensuring that the results are consistent and reproducible. The validation process involves several steps, including comparing the energy surfaces, reference energies, structures, and magnetic moments.

The `validate_consistency.py` script first compares the energy surfaces generated by both implementations. It reads the energy data from the output files and calculates the differences between the energies for each magnetization direction. If the differences exceed a specified tolerance, the script flags the discrepancies and provides detailed information about the specific points where the differences occur. This step is crucial for identifying any issues with the energy calculations and ensuring that the results are consistent.

Next, the script compares the reference energies used in both implementations. The reference energy is typically the energy of the collinear calculation with the magnetization aligned along the z-axis. The script reads the reference energy from the output files and calculates the difference between the two implementations. If the difference is significant, it may indicate issues with the reference calculation or the energy convergence.

The script also compares the structures used in both implementations. It reads the POSCAR files and checks for any differences in the atomic positions, lattice parameters, and chemical formulas. This step ensures that the calculations are performed on the same atomic configuration, which is essential for accurate comparisons.

Finally, the script compares the magnetic moments assigned to the atoms. It reads the magnetic moment data from the output files and calculates the differences between the two implementations. If the differences are significant, it may indicate issues with the magnetic moment initialization or the non-collinear calculations.

The `validate_consistency.py` script also includes a comprehensive reporting feature that generates a detailed report of the validation results. The report includes information about the energy differences, reference energy differences, structure differences, and magnetic moment differences. This report is useful for identifying any issues and ensuring that the results are consistent and reliable.

**Section sources**
- [validate_consistency.py](file://4_mae/validate_consistency.py#L1-L358)
- [MAE.py](file://4_mae/MAE.py#L1-L1086)

## Input Configuration and Examples
The input configuration for the MAE calculation module is defined in the `input_example_enhanced.py` file, which provides a comprehensive example of all available options for ensuring consistency between the original `MAE.py` script and the ported modular implementation. This file includes detailed documentation and usage instructions, making it easy for users to configure the calculations according to their specific needs.

The `input_example_enhanced.py` file defines several key parameters, including the number of grid points for the angles θ and φ, the VASP parameters, and the magnetic moment source. The number of grid points for each angle is configurable, allowing users to balance computational efficiency and accuracy. The VASP parameters, such as the energy cutoff (ENCUT), k-point mesh, and other settings, are also configurable, ensuring that the calculations are performed with the desired level of accuracy.

One of the key features of the `input_example_enhanced.py` file is the ability to specify the magnetic moment source. Users can choose between using element-based magnetic moments, which are based on the element type, or calculation-based magnetic moments, which are derived from the results of previous collinear calculations. This flexibility allows users to ensure that the magnetic moments are consistent with the original `MAE.py` script or the ported modular implementation.

The file also includes options for handling the structure, such as using the original structure or the primitive structure. This is important for ensuring that the calculations are performed on the correct atomic configuration. Additionally, the file includes options for using k-point optimized reference energy calculations, which can improve the accuracy of the results.

The `input_example_enhanced.py` file provides detailed usage instructions and troubleshooting guidelines, making it easy for users to configure the calculations and identify any issues. For example, users can enable the `use_mae_py_compatibility` option to use the exact same parameters as the original `MAE.py` script, ensuring that the results are consistent. Alternatively, users can selectively enable specific compatibility options to fine-tune the calculations.

**Section sources**
- [input_example_enhanced.py](file://4_mae/input_example_enhanced.py#L1-L192)
- [input_template.py](file://4_mae/input_template.py#L1-L124)

## Integration with FireWorks and VASP
The MAE calculation module in automag-1 is integrated with the FireWorks workflow management system and the VASP DFT code, ensuring efficient and reliable execution of the calculations. The integration with FireWorks allows for the automated submission and management of VASP jobs, making it easy to handle large-scale calculations on high-performance computing clusters. The `SubmitFirework.py` and `SubmitManual.py` scripts are used to submit the VASP jobs to the FireWorks database or to run them manually, depending on the user's preference.

The `SubmitFirework.py` script is responsible for submitting the VASP jobs to the FireWorks database. It creates a workflow that includes the necessary steps for the non-collinear DFT calculations, such as setting up the input files, running the VASP jobs, and processing the output. The script uses the `LaunchPad` class from the FireWorks library to connect to the remote database and submit the workflows. This integration ensures that the calculations are managed efficiently and that the results are stored in a centralized location.

The `SubmitManual.py` script, on the other hand, is used to submit the VASP jobs manually. This script is useful for users who prefer to run the calculations without the overhead of a workflow management system. The script creates the necessary input files and job scripts, and then submits the VASP jobs using the `sbatch` command. This approach provides more control over the job submission process and is suitable for smaller-scale calculations.

The integration with VASP is achieved through the `run_vasp.py` script, which is responsible for setting up the VASP environment and running the VASP executable. The script loads the necessary modules and libraries, such as the Intel MKL and MPI libraries, and then calls the `vasp_std` executable. The `run_vasp.py` script is configured to work with the specific VASP version and environment on the user's system, ensuring that the calculations are performed correctly.

The `run_vasp.py` script also includes options for specifying the VASP command and environment activation. Users can specify the command to run VASP, such as `mpirun vasp_std`, and the commands to activate and deactivate the Python environment. This flexibility allows users to customize the VASP execution to their specific needs.

The integration with FireWorks and VASP ensures that the MAE calculations are performed efficiently and reliably, making it easy to handle large-scale calculations on high-performance computing clusters. The combination of automated workflow management and manual job submission provides users with the flexibility to choose the approach that best suits their needs.

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L1-L260)
- [SubmitManual.py](file://common/SubmitManual.py#L1-L870)
- [run_vasp.py](file://ase/run_vasp.py#L1-L22)

## Troubleshooting Common Issues
Common issues in the MAE calculation module often arise from convergence failures in non-collinear calculations and inconsistent results. Convergence failures can occur due to various reasons, such as insufficient k-point mesh, inadequate energy cutoff (ENCUT), or issues with the magnetic moment initialization. To address these issues, users should first ensure that the k-point mesh and ENCUT values are sufficiently high to achieve convergence. The `1_submit.py` and `2_plot_results.py` scripts can be used to perform convergence tests and identify the optimal parameters.

Inconsistent results can also be caused by differences in the magnetic moment initialization. Users should verify that the magnetic moments are consistent with the original `MAE.py` script or the ported modular implementation. The `validate_consistency.py` script can be used to compare the magnetic moments and identify any discrepancies. If the magnetic moments are inconsistent, users should check the input configuration and ensure that the correct magnetic moment source is being used.

Another common issue is the failure of the VASP jobs to complete successfully. This can be due to issues with the job submission, such as incorrect job scripts or environment settings. Users should check the job scripts and ensure that the necessary modules and libraries are loaded correctly. The `run_vasp.py` script can be used to verify the VASP environment and ensure that the VASP executable is called correctly.

If the VASP jobs are failing due to memory or computational resource issues, users should consider increasing the number of CPU cores or the amount of memory allocated to the jobs. The `jobheader` parameter in the input configuration can be used to specify the number of CPU cores and the memory requirements. Additionally, users can use the `parallel_over_configurations` parameter to run the calculations in parallel, which can improve the efficiency of the calculations.

Finally, users should ensure that the input files and job scripts are correctly formatted and that all necessary files are present. The `input_example_enhanced.py` file provides a comprehensive example of the input configuration, which can be used as a reference for setting up the calculations. If issues persist, users can consult the `IMPLEMENTATION_SUMMARY.md` file for additional guidance and troubleshooting tips.

**Section sources**
- [MAE.py](file://4_mae/MAE.py#L1-L1086)
- [validate_consistency.py](file://4_mae/validate_consistency.py#L1-L358)
- [input_example_enhanced.py](file://4_mae/input_example_enhanced.py#L1-L192)
- [IMPLEMENTATION_SUMMARY.md](file://4_mae/IMPLEMENTATION_SUMMARY.md#L1-L223)