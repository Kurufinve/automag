# FireWorks-Based Orchestration

<cite>
**Referenced Files in This Document**   
- [SubmitFirework.py](file://common/SubmitFirework.py)
- [utilities.py](file://common/utilities.py)
- [1_submit.py](file://0_conv_tests/1_submit.py)
- [1_submit.py](file://1_lin_response/1_submit.py)
- [1_submit.py](file://4_mae/1_submit.py)
</cite>

## Table of Contents
1. [Introduction](#introduction)
2. [SubmitFirework Class Implementation](#submitfirework-class-implementation)
3. [Workflow Construction and Dependencies](#workflow-construction-and-dependencies)
4. [Firetasks in the Workflow Pipeline](#firetasks-in-the-workflow-pipeline)
5. [Convergence Testing Workflows](#convergence-testing-workflows)
6. [Perturbation Calculations](#perturbation-calculations)
7. [Links Dictionary and Job Dependencies](#links-dictionary-and-job-dependencies)
8. [Common Issues and Error Handling](#common-issues-and-error-handling)
9. [Performance Considerations](#performance-considerations)
10. [FireWorks Integration and Configuration](#fireworks-integration-and-configuration)

## Introduction
The FireWorks-based orchestration system provides a robust framework for managing computational workflows in materials science simulations. This documentation details the implementation of the SubmitFirework class, which serves as the primary interface for submitting FireWorks to a remote database. The system enables automated execution of VASP calculations with proper dependency management and state tracking. The orchestration framework supports various calculation modes including convergence testing, perturbation analysis, and single-point energy calculations. By leveraging the FireWorks workflow engine, the system ensures reliable execution of complex computational pipelines while providing mechanisms for error recovery and workflow debugging.

## SubmitFirework Class Implementation

The SubmitFirework class is designed to handle the submission of FireWorks to a remote database, providing a structured approach to workflow management. The class initialization accepts several key parameters that define the computational workflow:

- **poscar_file**: Specifies the path to the POSCAR file containing the crystal structure
- **mode**: Determines the calculation type (encut, kgrid, perturbations, singlepoint)
- **fix_params**: Dictionary containing fixed VASP calculation parameters
- **magmoms**: List of initial magnetic moments for the calculation
- **encut_values**: Range of ENCUT values for convergence testing
- **sigma_values**: Range of SIGMA values for k-point convergence
- **kpts_values**: Range of KPOINTS values for convergence testing
- **pert_values**: Values for perturbation calculations
- **name**: Custom name for the workflow
- **dummy_atom**: Atom type to be used as dummy in perturbation calculations
- **dummy_position**: Position index of the dummy atom in the structure

The class constructor validates parameter combinations based on the specified mode, ensuring that only appropriate parameter sets are used for each calculation type. For example, when mode is set to 'encut', the constructor verifies that encut_values is provided while other variable parameter sets are None. The initialization also sets up the LaunchPad connection using a YAML configuration file, establishing the database connection for workflow submission.

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L26-L77)

## Workflow Construction and Dependencies

The workflow construction process in the SubmitFirework class follows a systematic approach to create interconnected FireWorks with proper dependencies. The add_wflow method is responsible for building the workflow structure based on the specified calculation mode. For most modes, the workflow begins with a single-point calculation followed by a recalculation step that uses magnetic moments from the previous run.

The workflow construction process involves several key steps:
1. Reading the atomic structure from the POSCAR file
2. Creating encoded representation of the structure
3. Building Firework objects for each calculation step
4. Establishing dependencies between Fireworks using the links_dict
5. Submitting the complete workflow to the LaunchPad

For perturbation calculations, the workflow includes additional steps: non-self-consistent (NSC) and self-consistent (SC) calculations for each perturbation value, followed by charge writing tasks. The workflow structure adapts dynamically based on the calculation mode, ensuring that appropriate steps are included for each type of analysis.

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L104-L258)

## Firetasks in the Workflow Pipeline

The workflow pipeline utilizes specialized Firetasks to execute specific computational steps. The primary Firetasks include:

### VaspCalculationTask
This Firetask handles the execution of VASP calculations with configurable parameters. It supports various calculation modes including magnetic calculations and perturbation analysis. The task can accept an encoded structure or read the input from the previous step's output. For magnetic calculations, it properly sets up the ISIPN and LORBIT parameters and initializes magnetic moments. The task also handles convergence monitoring by writing convergence status to a file.

### WriteOutputTask
Responsible for collecting and formatting calculation results, this Firetask generates output files containing energy, convergence status, chemical formula, space group, and magnetic moments. It reads information from multiple sources including OUTCAR, CONTCAR, and job specification. The task supports both energy and enthalpy output with appropriate corrections for convergence errors.

### WriteChargesTask
Specifically designed for perturbation calculations, this Firetask extracts and writes charge information from VASP output files. It identifies the dummy atom position and extracts f-electron or d-electron charges based on the calculation type. The task includes validation to ensure magnetic moments remain within acceptable ranges before writing charge data.

**Section sources**
- [utilities.py](file://common/utilities.py#L64-L151)
- [utilities.py](file://common/utilities.py#L155-L218)

## Convergence Testing Workflows

The system supports two primary types of convergence testing: energy cutoff (ENCUT) and k-point grid (kgrid) convergence. These workflows are demonstrated in the 0_conv_tests module, where the 1_submit.py script shows practical implementation examples.

For ENCUT convergence testing, the workflow iterates through a range of energy cutoff values, typically from 500 to 1000 eV in 10 eV increments. Each calculation produces a single-point energy that is used to assess convergence behavior. The workflow automatically generates appropriate names for each Firework based on the ENCUT value.

For k-point convergence testing, the system evaluates both SIGMA and KPOINTS parameters simultaneously. The workflow creates calculations for various combinations of these parameters, allowing comprehensive assessment of k-point convergence. The naming convention combines both parameter values to uniquely identify each calculation.

The convergence testing workflows include automatic output generation that collects results in a structured format, facilitating subsequent analysis and plotting. The system also handles magnetic configuration determination, automatically assigning appropriate magnetic moments based on element types when not explicitly specified.

**Section sources**
- [1_submit.py](file://0_conv_tests/1_submit.py#L79-L97)

## Perturbation Calculations

Perturbation calculations are implemented through the 'perturbations' mode in the SubmitFirework class. This mode enables systematic analysis of electronic structure changes by introducing a dummy atom with varying perturbation values. The workflow for perturbation calculations follows a multi-step process:

1. **Structure Modification**: The specified atom position is replaced with a dummy atom (e.g., Zn) to create a perturbed system
2. **Base Calculation**: A single-point calculation establishes the reference state
3. **NSC Calculations**: Non-self-consistent calculations for each perturbation value, using WAVECAR and CHGCAR from the base calculation
4. **SC Calculations**: Self-consistent calculations for each perturbation value
5. **Charge Extraction**: Writing charge differences for analysis

The implementation in 1_lin_response/1_submit.py demonstrates this workflow, specifying parameters such as the dummy atom type, position, and perturbation values. The system automatically handles the creation of appropriate VASP input parameters, including LDAU settings for the dummy atom. The workflow ensures proper dependency ordering, with NSC calculations depending on the base calculation and SC calculations depending on their corresponding NSC results.

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L150-L210)
- [1_submit.py](file://1_lin_response/1_submit.py#L60-L77)

## Links Dictionary and Job Dependencies

The links_dict mechanism is central to establishing proper job dependencies and execution order in the FireWorks workflow. This dictionary maps Firework IDs to their dependent Fireworks, ensuring sequential execution according to the workflow logic.

The dependency structure varies based on the calculation mode:
- **Standard workflows**: Single-point calculation (fw_id=0) → Recalculation (fw_id=1) → Output writing (fw_id=2)
- **Perturbation workflows**: Base calculation → Multiple NSC calculations → Multiple SC calculations → Multiple charge writing tasks

The links_dict is constructed by analyzing the fireworks list structure, where each level represents a stage in the workflow. The algorithm handles different connection patterns:
- One-to-many: Single Firework depending on multiple subsequent Fireworks
- Many-to-one: Multiple Fireworks depending on a single subsequent Firework
- Many-to-many: Corresponding Fireworks at each level depending on each other

This flexible dependency system ensures that calculations proceed in the correct order while allowing parallel execution of independent tasks where appropriate. The implementation automatically generates unique Firework IDs and establishes the appropriate links between them.

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L220-L250)

## Common Issues and Error Handling

The system addresses several common issues that may arise during workflow execution:

### Workflow Submission Failures
Submission failures can occur due to improper parameter combinations or missing required parameters. The SubmitFirework constructor includes comprehensive validation to catch these issues early. For example, attempting to run an 'encut' mode calculation without specifying encut_values will raise a ValueError with a descriptive message.

### Database Connection Problems
The LaunchPad connection is established at the module level using a YAML configuration file. Connection issues typically stem from incorrect file paths or database authentication problems. The system uses a hardcoded path to the launchpad file, which should be updated to match the user's environment.

### Handling Failed Fireworks
The system includes mechanisms to identify and handle failed calculations:
- Convergence status is written to 'is_converged' files for each calculation
- The WriteOutputTask reads these files to report convergence status
- Failed calculations can be identified by 'NONCONVERGED' entries
- Manual intervention may be required to restart failed calculations

Error recovery strategies include:
- Verifying input parameters before submission
- Checking for existing completed calculations to avoid duplication
- Using appropriate VASP parameters to improve convergence
- Implementing proper resource allocation in job scripts

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L252-L258)
- [utilities.py](file://common/utilities.py#L145-L151)

## Performance Considerations

When submitting large-scale workflows, several performance considerations should be addressed:

### Resource Management
Large workflows can consume significant computational resources. The system should be configured to:
- Limit concurrent job submissions to avoid overwhelming the cluster
- Use appropriate walltime and resource allocations in job scripts
- Implement checkpointing for long-running calculations

### Workflow Optimization
To improve efficiency:
- Use appropriate k-point sampling to balance accuracy and computational cost
- Optimize ENCUT values based on preliminary tests
- Consider using symmetry operations to reduce computational load
- Implement intelligent workflow branching based on intermediate results

### Submission Strategies
For large parameter spaces:
- Submit calculations in batches rather than all at once
- Use priority queues to manage job scheduling
- Implement monitoring to detect and address stalled calculations
- Consider using workflow-level parallelization where appropriate

The system's modular design allows for easy adaptation to different performance requirements and computational environments.

## FireWorks Integration and Configuration

The integration with the FireWorks LaunchPad system requires proper configuration of the database connection. The current implementation uses a hardcoded path to the launchpad YAML file:

```python
launchpad = LaunchPad.from_file('/home/mgalasso/.fireworks/my_launchpad.yaml')
```

This configuration should be updated to match the user's environment. The launchpad file typically contains database connection parameters including:
- Hostname and port
- Database name
- Username and authentication credentials
- SSL configuration

The system assumes that the FireWorks database is properly set up and accessible from the submission environment. Users should verify that:
- The MongoDB server is running and accessible
- The database user has appropriate permissions
- Network connectivity is established between the client and server
- Firewall rules allow the necessary connections

Remote database submission requires that the FireWorks rocket executable is available on the target compute nodes and that the database connection parameters are correctly configured in the launchpad file.

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L24-L25)