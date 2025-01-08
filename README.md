# OptimaNBS

## Table Of Contents

- [Description](#description)
- [Java Code Usage](#java-code-usage)
- [Python Code Usage](#python-code-usage)
- [Contribution](#contribution)
- [Citation](#citation)
- [License](#license)

## Description
Considering the relevance of Operational Research (OR) in urban planning, this project releases the source code implemented with reference to the manuscript "Optimal Placement of Nature-Based Solutions for Urban Challenges" to exploit OR tools and develop a flexible mathematical formulation for locating multiple Nature-Based Solutions (NBSs) while considering real-world constraints. Specifically, this code was used for the following contributions:

- **Implementation of a Mixed-Integer Linear Programming (MILP) model:** This model optimizes the placement of various NBSs within urban environments, addressing metrics of different Urban Challenges (UCs), including ground temperature reduction, air quality improvement, and enhanced accessibility to green areas.
- **Evaluation of the MILP formulation:** The model's effectiveness is tested using real-world case studies from Italian cities, selected for their diverse urban characteristics, climatic conditions, and geographic contexts.
- **Instances creation and Impact Estimation:** Data on both NBSs and UCs are gathered from state-of-the-art literature and integrated with national databases (e.g., population and land use) to create a set of benchmark instances for comprehensive testing.

## Java Code Usage
Follow these steps to use the Java code:

1. Download the content of the `Code/Java` folder.
2. Open it with an appropriate IDE (IntelliJ IDEA was used during development).
3. Integrate information related to the Java version you intend to use.
4. Link the available Gurobi license and add a reference to the Gurobi JAR file (check the Gurobi website for installation and usage instructions).
5. Run the code contained in the `Main.java` class.

### Additional Notes
- The `Code/Java/Test` folder contains an example instance for testing the method. The results will be saved in the `Code/Java/Results` folder.
- The `Code/Java/Config` folder contains the kernel used during the testing phase.
- All parameters can be customized by accessing and modifying the content of `Code/Java/src/main/java/Config.java`.

## Python Code Usage
Follow these steps to use the Python code:

1. Download the content of the `Code/Python` folder.
2. Open it with an appropriate IDE (PyCharm was used during development).
3. Install the required dependencies.
4. To create the kernel files, run the `Code/Python/Kernels/KernelUtils.py` script.
5. To generate instances for model testing, run the `Code/Python/InstanceGenerator/InstanceGeneratorUtils.py` script.

### Additional Notes
- Both procedures can be customized by modifying the parameters in the files `Code/Python/InstanceGenerator/config_instance_generator.py` and `Code/Python/Kernels/config_kernel_generator.py`.
- The `Code/Python/Data` folder contains the necessary input data for instance generation, and output data will be stored in the same folder after execution.

## Contribution
This project is the result of research conducted by the following three researchers:

1. [Diego Maria Pinto - diegomaria.pinto@iasi.cnr.it - Institute for Systems Analysis and Computer Science "Antonio Ruberti", National Research Council, Via dei Taurini 19, 00185 Rome, Italy]
2. [Davide Donato Russo - davide.donatorusso@iasi.cnr.it - Institute for Systems Analysis and Computer Science "Antonio Ruberti", National Research Council, Via dei Taurini 19, 00185 Rome, Italy]
3. [Antonio M. Sudoso - antoniomaria.sudoso@uniroma1.it - Department of Computer, Control and Management Engineering "Antonio Ruberti", Sapienza University of Rome, Via Ariosto 25, 00185 Rome, Italy - Institute for Systems Analysis and Computer Science "Antonio Ruberti", National Research Council, Via dei Taurini 19, 00185 Rome, Italy]

Feel free to contact us for potential collaborations.

## Citation
You are free to use this work, but remember to cite the correct reference. To cite this work, please use the following reference:

```
[TODO: INSERT REFERENCE]
```

## License
This project is licensed under the GNU General Public License (GPL). See the `LICENSE` file for more details.

