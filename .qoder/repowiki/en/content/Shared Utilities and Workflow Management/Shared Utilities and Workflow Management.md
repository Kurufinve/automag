# Shared Utilities and Workflow Management

<cite>
**Referenced Files in This Document**   
- [SubmitFirework.py](file://common/SubmitFirework.py)
- [SubmitManual.py](file://common/SubmitManual.py)
- [utilities.py](file://common/utilities.py)
- [write_charges.py](file://common/write_charges.py)
- [write_output.py](file://common/write_output.py)
- [write_magmoms_input.py](file://common/write_magmoms_input.py)
- [get_magmoms_vasp.py](file://common/get_magmoms_vasp.py)
- [write_output_vasp.py](file://common/write_output_vasp.py)
</cite>

## Table of Contents
1. [Introduction](#introduction)
2. [Workflow Orchestration](#workflow-orchestration)
3. [Data Processing Pipelines](#data-processing-pipelines)
4. [Output Formatting Standards](#output-formatting-standards)
5. [Component Interactions](#component-interactions)
6. [Hybrid Execution Models](#hybrid-execution-models)
7. [Infrastructure Requirements](#infrastructure-requirements)
8. [System Context Diagrams](#system-context-diagrams)
9. [Cross-Cutting Concerns](#cross-cutting-concerns)

## Introduction
The shared utilities component provides a comprehensive framework for managing computational materials science workflows, particularly focused on magnetic property calculations. This system implements dual execution paradigms through FireWorks-based automated workflow management and manual job submission modes. The architecture enables systematic exploration of material properties through convergence testing, perturbation analysis, and magnetic anisotropy calculations. The component integrates tightly with VASP (Vienna Ab initio Simulation Package) for electronic structure calculations while providing robust data processing, error handling, and result aggregation capabilities.

## Workflow Orchestration
The workflow management system provides two complementary approaches for executing computational workflows: automated FireWorks orchestration and manual job submission. The FireWorks-based approach enables distributed, fault-tolerant execution of complex calculation sequences, while the manual mode provides direct control over job submission for environments without workflow management infrastructure.

```mermaid
graph TD
A[User Input] --> B{Execution Mode}
B --> C[FireWorks Mode]
B --> D[Manual Mode]
C --> E[SubmitFirework]
D --> F[SubmitManual]
E --> G[FireWorks Database]
F --> H[Job Scripts]
G --> I[Remote Execution]
H --> J[Local/Cluster Execution]
I --> K[Result Collection]
J --> K
K --> L[Data Processing]
L --> M[Output Generation]
```

**Diagram sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L25-L258)
- [SubmitManual.py](file://common/SubmitManual.py#L24-L853)

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L25-L258)
- [SubmitManual.py](file://common/SubmitManual.py#L24-L853)

## Data Processing Pipelines
The system implements specialized data processing pipelines for magnetic moments and charge analysis, with different implementations for FireWorks and manual execution modes. These pipelines extract, validate, and transform raw calculation outputs into scientifically meaningful results.

### Magnetic Moments Processing
The magnetic moments pipeline follows a two-stage approach: initial moment extraction from converged calculations and subsequent validation against convergence criteria. The system captures both initial and final magnetic moments, enabling analysis of magnetic configuration evolution during self-consistent field calculations.

```mermaid
flowchart TD
A[VASP Calculation] --> B{Convergence Check}
B --> |Converged| C[Extract Magnetic Moments]
B --> |Not Converged| D[Flag for Review]
C --> E[Validate Moment Stability]
E --> F[Store Initial/Final Moments]
F --> G[Generate Output File]
```

**Diagram sources**
- [get_magmoms_vasp.py](file://common/get_magmoms_vasp.py#L1-L42)
- [write_magmoms_input.py](file://common/write_magmoms_input.py#L1-L50)

### Charges Processing Pipeline
The charge analysis pipeline focuses on extracting charge information for perturbation calculations, particularly for Hubbard U parameter determination. The system implements validation logic to ensure result reliability by checking magnetic moment stability and convergence status.

```mermaid
flowchart TD
A[NSC Calculation] --> B[Extract Charges]
A --> C[Check Convergence]
A --> D[Compare Magnetic Moments]
B --> E{Valid Results?}
C --> E
D --> E
E --> |Yes| F[Write Charges]
E --> |No| G[Discard Results]
```

**Diagram sources**
- [write_charges.py](file://common/write_charges.py#L1-L125)
- [utilities.py](file://common/utilities.py#L221-L276)

**Section sources**
- [write_charges.py](file://common/write_charges.py#L1-L125)
- [utilities.py](file://common/utilities.py#L221-L276)

## Output Formatting Standards
The system implements consistent output formatting standards across different calculation types and execution modes. Output files contain structured information including calculation metadata, convergence status, chemical formula, crystallographic information, and physical properties.

### Output Structure
Output files follow a standardized format with the following fields:
- **System identifier**: Unique identifier for the calculation
- **Firework ID**: Workflow tracking identifier (FireWorks mode)
- **Convergence status**: Per-step convergence information
- **Chemical formula**: Standardized chemical formula
- **Space group**: Crystallographic space group symbol
- **Energy/Enthalpy**: Total energy or enthalpy with correction
- **Magnetic moments**: Initial and final magnetic moments (when applicable)

```mermaid
erDiagram
OUTPUT_FILE {
string system_id PK
int firework_id
string convergence_status
string chemical_formula
string space_group
float energy
float enthalpy
string initial_magmoms
string final_magmoms
}
```

**Diagram sources**
- [write_output.py](file://common/write_output.py#L1-L147)
- [write_output_vasp.py](file://common/write_output_vasp.py#L1-L113)

**Section sources**
- [write_output.py](file://common/write_output.py#L1-L147)
- [write_output_vasp.py](file://common/write_output_vasp.py#L1-L113)

## Component Interactions
The system architecture features well-defined interactions between core components, enabling modular design and separation of concerns. The SubmitFirework and VaspCalculationTask components form the backbone of the automated workflow system, while supporting utilities handle data processing and output generation.

### SubmitFirework and VaspCalculationTask Interaction
The SubmitFirework class orchestrates workflow creation by instantiating VaspCalculationTask instances for each calculation step. This interaction enables parameterized execution of VASP calculations within the FireWorks framework.

```mermaid
classDiagram
class SubmitFirework {
+__init__(poscar_file, mode, fix_params, magmoms, ...)
+submit()
+add_wflow(params, name)
}
class VaspCalculationTask {
+run_task(fw_spec)
}
class WriteOutputTask {
+run_task(fw_spec)
}
class WriteChargesTask {
+run_task(fw_spec)
}
SubmitFirework --> VaspCalculationTask : "creates"
SubmitFirework --> WriteOutputTask : "creates"
SubmitFirework --> WriteChargesTask : "creates"
VaspCalculationTask --> Firework : "extends"
WriteOutputTask --> Firework : "extends"
WriteChargesTask --> Firework : "extends"
```

**Diagram sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L25-L258)
- [utilities.py](file://common/utilities.py#L64-L276)

### Structure Encoding Utilities
The system implements utilities for encoding and decoding atomic structures to facilitate workflow persistence and data exchange. The atoms_to_encode and encode_to_atoms functions provide bidirectional conversion between ASE Atoms objects and JSON-serializable representations.

```mermaid
classDiagram
class utilities {
+atoms_to_encode(atoms)
+encode_to_atoms(encode)
}
class Atoms {
+get_cell()
+get_scaled_positions()
+get_atomic_numbers()
+get_pbc()
}
utilities --> Atoms : "reads from"
Atoms --> utilities : "provides data"
```

**Diagram sources**
- [utilities.py](file://common/utilities.py#L1-L62)

**Section sources**
- [utilities.py](file://common/utilities.py#L1-L62)

## Hybrid Execution Models
The system implements a hybrid execution model that supports both automated workflow management through FireWorks and manual job submission. This dual approach provides flexibility for different computational environments and user requirements.

### Technical Implementation
The hybrid model is implemented through parallel class hierarchies: SubmitFirework for automated execution and SubmitManual for script-based submission. Both classes share similar interfaces and configuration parameters, enabling consistent workflow definition across execution modes.

```mermaid
classDiagram
class SubmitFirework {
+submit()
+add_wflow()
}
class SubmitManual {
+submit()
+write_input()
+write_jobscript()
}
SubmitFirework --> LaunchPad : "submits to"
SubmitManual --> ShellScript : "generates"
SubmitFirework --> Firework : "creates"
SubmitManual --> Vasp : "configures"
```

**Diagram sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L25-L258)
- [SubmitManual.py](file://common/SubmitManual.py#L24-L853)

### Error Handling in Job Submission
The system implements comprehensive error handling for job submission, with different strategies for the two execution modes. The FireWorks mode leverages the framework's built-in fault tolerance, while the manual mode implements validation checks and recovery procedures.

```mermaid
flowchart TD
A[Job Submission] --> B{Execution Mode}
B --> C[FireWorks]
B --> D[Manual]
C --> E[Queue Management]
C --> F[Automatic Retry]
C --> G[Error Logging]
D --> H[Script Validation]
D --> I[Directory Creation]
D --> J[File Existence Check]
E --> K[Execution]
F --> K
G --> K
H --> K
I --> K
J --> K
K --> L[Result Monitoring]
```

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L25-L258)
- [SubmitManual.py](file://common/SubmitManual.py#L24-L853)

## Infrastructure Requirements
The system has specific infrastructure requirements that differ between the two execution modes, particularly regarding database and queue management systems.

### FireWorks Mode Requirements
The automated workflow mode requires:
- **MongoDB**: Database backend for FireWorks workflow persistence
- **FireWorks server**: Workflow management service
- **Queue adapter**: Integration with HPC job schedulers (SLURM, PBS, etc.)
- **LaunchPad configuration**: YAML file with database connection parameters

### Manual Mode Requirements
The script-based execution mode requires:
- **Shell environment**: Bash or compatible shell
- **Job scheduler access**: SLURM, PBS, or similar (optional)
- **Directory structure**: Properly configured calculation directories
- **Environment modules**: Python and VASP dependencies

```mermaid
graph TD
A[Execution Mode] --> B[FireWorks]
A --> C[Manual]
B --> D[MongoDB]
B --> E[FireWorks Server]
B --> F[Queue Adapter]
B --> G[LaunchPad]
C --> H[Shell Environment]
C --> I[Job Scheduler]
C --> J[Directory Structure]
C --> K[Environment Modules]
```

**Section sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L1-L25)
- [SubmitManual.py](file://common/SubmitManual.py#L1-L23)

## System Context Diagrams
The system context diagrams illustrate the complete workflow generation, execution, and result collection process for both execution modes.

### FireWorks Workflow Context
```mermaid
sequenceDiagram
participant User
participant SubmitFirework
participant FireWorksDB
participant HPC
participant ResultProcessor
User->>SubmitFirework : Configure workflow
SubmitFirework->>FireWorksDB : Submit workflow
FireWorksDB->>HPC : Execute jobs
HPC-->>FireWorksDB : Report completion
FireWorksDB->>ResultProcessor : Trigger output generation
ResultProcessor-->>User : Deliver results
```

**Diagram sources**
- [SubmitFirework.py](file://common/SubmitFirework.py#L25-L258)

### Manual Execution Context
```mermaid
sequenceDiagram
participant User
participant SubmitManual
participant Script
participant HPC
participant ResultProcessor
User->>SubmitManual : Configure calculation
SubmitManual->>Script : Generate job script
Script->>HPC : Submit job
HPC-->>Script : Execute calculation
Script->>ResultProcessor : Process output
ResultProcessor-->>User : Deliver results
```

**Diagram sources**
- [SubmitManual.py](file://common/SubmitManual.py#L24-L853)

## Cross-Cutting Concerns
The system addresses several cross-cutting concerns that impact multiple components and execution modes.

### Logging and Monitoring
The system implements comprehensive logging at multiple levels:
- **Workflow logging**: FireWorks job information and execution status
- **Calculation logging**: VASP output and convergence information
- **Error logging**: Exception handling and failure reporting
- **Progress logging**: Step-by-step execution tracking

### Error Recovery
The system implements several error recovery mechanisms:
- **Automatic retry**: FireWorks workflow resilience
- **Checkpointing**: State preservation between calculation steps
- **Validation checks**: Pre-execution verification of inputs and directories
- **Fallback procedures**: Alternative execution paths when primary methods fail

### HPC Environment Compatibility
The system is designed to work across different HPC environments through:
- **Configurable job headers**: Support for SLURM, PBS, and other schedulers
- **Environment activation/deactivation**: Module loading and unloading
- **Path configuration**: Environment variable-based path resolution
- **Flexible directory structure**: Adaptable to different filesystem layouts

```mermaid
graph TD
A[HPC Environment] --> B[Job Scheduler]
A --> C[Module System]
A --> D[Filesystem]
B --> E[SLURM]
B --> F[PBS]
B --> G[LSF]
C --> H[Environment Modules]
C --> I[Lmod]
D --> J[NFS]
D --> K[Lustre]
E --> SubmitManual
F --> SubmitManual
G --> SubmitManual
H --> SubmitManual
I --> SubmitManual
J --> System
K --> System
```

**Section sources**
- [SubmitManual.py](file://common/SubmitManual.py#L24-L853)