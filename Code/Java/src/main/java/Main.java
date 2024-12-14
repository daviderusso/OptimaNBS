import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

public class Main {
    public static void main(String[] args) {
        // Initialize configuration object
        Config config = new Config();

        // Create a File object pointing to the folder with test instances
        File folder = new File(Config.TestInstanceFolder);
        // Get all files in the folder
        File[] listOfFiles = folder.listFiles();

        // Iterate through each file in the folder
        for (File file : listOfFiles) {
            if (file.isFile()) { // Process only files (not directories)
                try {
                    // Output the file name
                    System.out.println(file.getName() + "##################################");

                    // Read Kernel configuration file
                    System.out.println("Read Kernel ##################################");
                    String pathKernel = Config.KernelFolder + Config.KernelName;
                    // Read kernels and store in a nested HashMap
                    HashMap<String, HashMap<String, Kernel>> kernelSet = KernelUtils.readKernels(pathKernel);

                    // Read the specific instance from the current file
                    System.out.println("Read Instance ##################################");
                    Instance instance = InstanceUtils.readInstance(file.getPath());
                    // Compute the maximum budget for the current instance based on its width and height
                    config.budgetComputationMax(instance.W, instance.H);

                    //--------------------------------------- INIT FF AND OUT
                    // Compute the initial Objective Function (OF) values for the instance
                    HashMap<String, Double> initFF = InstanceUtils.computeOF(instance, config, instance.A, new int[instance.W][instance.H], kernelSet);
                    // Output initial information about the results
                    ResultsUtils.getInfos(file.getName(), instance, Config.OutputFolder, initFF);

                    //--------------------------------------- Normalize Data matrix
                    System.out.println("Normalize Data ##################################");
                    // Normalize the instance data and kernel
                    instance = InstanceUtils.normalizeData(instance, config);
                    kernelSet = KernelUtils.normalizeKernel(instance, kernelSet);

                    //--------------------------------------- CLUSTERING
                    System.out.println("Clusterize Park Space ##################################");
                    // Perform clustering on park spaces in the instance
                    int[][] clusteredMapPark = InstanceUtils.getParkClusterCCLWithBounds(instance, config);
                    // Compute the size of each cluster
                    ArrayList<Integer> sizeCluster = InstanceUtils.computeClusterSize(clusteredMapPark);

                    //--------------------------------------- INITIAL FF AND REDUCTION
                    System.out.println("Compute Initial FF ##################################");
                    // Compute the maximum value for the objective function without applying Theta
                    HashMap<String, Double> singleFF = InstanceUtils.computeMaxFFNoTheta(instance, instance.A, new int[instance.W][instance.H], kernelSet);
                    // Compute fairness bounds for urban challenges
                    double[] boundFairness = InstanceUtils.computeBoundFairness(instance, kernelSet.get(instance.UrbanChallenges[instance.UrbanChallenges.length - 1]));
                    // Arrays to hold the max values, reference values, and reduction percentages
                    double[] maxVal = new double[singleFF.size()]; // MAX VALUE OF THE MATRIX
                    double[] refVal = new double[singleFF.size()]; // MAX REDUCABLE VALUE IN THE MATRIX
                    double[] deltaB = new double[singleFF.size()]; // REDUCTION PERCENTAGE
                    // For each urban challenge, set the values for max, reference, and reduction percentage
                    for (int b = 0; b < instance.UrbanChallenges.length; b++) {
                        maxVal[b] = singleFF.get(instance.UrbanChallenges[b]);
                        refVal[b] = singleFF.get(instance.UrbanChallenges[b]) * (1.0 - Config.percBoundUCs);
                        deltaB[b] = singleFF.get(instance.UrbanChallenges[b]) * Config.percBoundUCs;
                    }

                    // Compute the MILP model
                    System.out.println("Compute Model ##################################");
                    // Initialize and solve the MILP model using the Gurobi solver
                    SolverMILPGurobi solver = new SolverMILPGurobi(instance, kernelSet, Config.showOut, Config.timeLimit, refVal, deltaB, maxVal, clusteredMapPark, sizeCluster, boundFairness, config);
                    // Execute the MILP computation
                    solver.compute_MILP();
                    System.out.println("####################################################################");
                } catch (Exception e) {
                    // Handle any errors during processing and print the error details
                    System.err.println("ERROR INSTANCE: " + file.getName());
                    e.printStackTrace();
                }
            }
        }
    }
}
